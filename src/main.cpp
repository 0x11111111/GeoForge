#include <memory>
#include <sstream>

#include <vtkSelection.h>
#include <vtkInformation.h>
#include <vtkSelectionNode.h>
#include <vtkCellArray.h>
#include <vtkLineSource.h>
#include <vtkLookupTable.h>
#include <vtkFloatArray.h>

#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/AABB_tree.h>
#include <CGAL/AABB_traits.h>
#include <CGAL/AABB_face_graph_triangle_primitive.h>

//#include <igl/avg_edge_length.h>
//#include <igl/cotmatrix.h>
//#include <igl/invert_diag.h>
//#include <igl/massmatrix.h>
//#include <igl/parula.h>
//#include <igl/per_corner_normals.h>
//#include <igl/per_face_normals.h>
//#include <igl/per_vertex_normals.h>
//#include <igl/principal_curvature.h>
//#include <igl/read_triangle_mesh.h>

#include "Render.h"
#include "VectorAlgorithm.h"
#include "Rotation.h"
#include "Timer.hpp"
vtkSmartPointer<vtkActor> PolyDataActor;
vtkSmartPointer<vtkActor> CutPolyDataActor;
vtkSmartPointer<vtkActor> TeethPolyDataActor;
vtkSmartPointer<vtkActor> BitePolyDataActor;

typedef Kernel::Segment_3 Segment_3;

typedef CGAL::AABB_face_graph_triangle_primitive<SurfaceMesh> Primitive;
typedef CGAL::AABB_traits<Kernel, Primitive> Traits;
typedef CGAL::AABB_tree<Traits> Tree;
typedef Tree::Primitive_id Primitive_id;
typedef boost::optional<Tree::Intersection_and_primitive_id<Kernel::Ray_3>::Type> Ray_intersection;
typedef boost::optional<Tree::Intersection_and_primitive_id<Segment_3>::Type> Segment_intersection;
typedef Kernel::Ray_3 Ray;

int SurfaceMeshToPolyData(const SurfaceMesh& mesh, vtkSmartPointer<vtkPolyData>& polyData)
{
    //polyData->Initialize();
    polyData = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
    VTKPoints->SetNumberOfPoints(mesh.vertices().size());

#pragma omp parallel for
    for (int i = 0; i < mesh.vertices().size(); i++)
    {
        const Kernel::Point_3& point = mesh.point(vertex_descriptor(i));
        VTKPoints->SetPoint(i, point.x(), point.y(), point.z());
    }

    vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfComponents(4);
    connectivity->SetNumberOfTuples(mesh.faces().size());

    //#pragma omp parallel for
    for (int j = 0; j < mesh.faces().size(); j++)
    {
        int ids[3] = { 0 };
        int index = 0;
        for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
        {
            //cout << index<<" ";
            //cout << mesh.target(h).idx() << endl;
            ids[index] = mesh.target(h).idx();
            index++;
        }
        if (index > 3)
        {
            std::cerr << "SurfaceMeshToPolyData Failed" << std::endl;
        }

        connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
    }

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(mesh.faces().size(), connectivity);

    polyData->SetPoints(VTKPoints);
    polyData->SetPolys(cells);
    //polyData->Modified();

    return 0;

}

int SurfaceMeshToPolyDataWithColor(
    const SurfaceMesh& mesh, 
    vtkSmartPointer<vtkPolyData>& polyData, 
    const std::set<vertex_descriptor>& intersected_vertex_set,
    vtkSmartPointer<vtkFloatArray>& distances_color_array
)
{
    //polyData->Initialize();
    polyData = vtkSmartPointer<vtkPolyData>::New();

    vtkSmartPointer<vtkPoints> VTKPoints = vtkSmartPointer<vtkPoints>::New();
    VTKPoints->SetNumberOfPoints(mesh.vertices().size());
    VLMap distance_property_map = teeth_sm.property_map<vertex_descriptor, double>("v:distance_from_bite").first;
//#pragma omp parallel for
    for (int i = 0; i < mesh.vertices().size(); i++)
    {
        vertex_descriptor vd = vertex_descriptor(i);
        const Kernel::Point_3& point = mesh.point(vd);
        VTKPoints->SetPoint(i, point.x(), point.y(), point.z());

        if (intersected_vertex_set.find(vd) != intersected_vertex_set.end())
        {
            distances_color_array->InsertNextValue(distance_property_map[vd]);
		}
		else
		{
			distances_color_array->InsertNextValue(0.0);
        }
    }

    vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
    connectivity->SetNumberOfComponents(4);
    connectivity->SetNumberOfTuples(mesh.faces().size());

    //#pragma omp parallel for
    for (int j = 0; j < mesh.faces().size(); j++)
    {
        int ids[3] = { 0 };
        int index = 0;
        for (const halfedge_descriptor& h : mesh.halfedges_around_face(mesh.halfedge(face_descriptor(j))))
        {
            //cout << index<<" ";
            //cout << mesh.target(h).idx() << endl;
            ids[index] = mesh.target(h).idx();
            index++;
        }
        if (index > 3)
        {
            std::cerr << "SurfaceMeshToPolyData Failed" << std::endl;
        }

        connectivity->SetTuple4(j, 3, ids[0], ids[1], ids[2]);
    }

    vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
    cells->SetCells(mesh.faces().size(), connectivity);

    polyData->SetPoints(VTKPoints);
    polyData->SetPolys(cells);

    polyData->GetPointData()->SetScalars(distances_color_array);
    //polyData->Modified();

    return 0;

}

int PolyDataToSurfaceMesh(vtkPolyData* polyData, SurfaceMesh& surfaceMesh)
{
    vtkPoints* points = polyData->GetPoints();
    int numberOfPoints = polyData->GetNumberOfPoints();

    if (0 == numberOfPoints || points == NULL)
        return -1;

    for (int i = 0; i < numberOfPoints; i++)
    {
        double point[3];
        points->GetPoint(i, point);
        surfaceMesh.add_vertex(Kernel::Point_3(point[0], point[1], point[2]));
    }

    vtkCellArray* triangles = polyData->GetPolys();
    triangles->InitTraversal();

    //vtkIdType npts = 3, * cell;
    vtkIdType npts;               // 存储单元的点数
    vtkIdType* cell;         // 存储单元的点索引数组
    while (triangles->GetNextCell(npts, cell))
    {
        surfaceMesh.add_face(CGAL::SM_Vertex_index(cell[0]), CGAL::SM_Vertex_index(cell[1]), CGAL::SM_Vertex_index(cell[2]));
    }

    return 0;
}

void LeftPress(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

void MouseMove(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

void LeftRelease(vtkObject* caller, long unsigned int eventId, void* clientData, void* callData)
{
}

vtkSmartPointer<vtkActor> CreateNormalLineActor(const Point_3& point, const Vector_3& normal, double length = 1.0, const double r = 1.0, const double g = 1.0, const double b = 1.0)
{
    Point_3 offset(normal.x() * length, normal.y() * length, normal.z() * length);  
    Point_3 endPoint(point.x() + offset.x(), point.y() + offset.y(), point.z() + offset.z());

    vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
    lineSource->SetPoint1(point.x(), point.y(), point.z());
    lineSource->SetPoint2(endPoint.x(), endPoint.y(), endPoint.z());

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(lineSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(1, 1, 1);  // Red color for example
    actor->GetProperty()->SetLineWidth(1);

    return actor;
}

vtkSmartPointer<vtkActor> CreateNormalLineActor(const Point_3& point_a, const Point_3& point_b, const double r = 1.0, const double g = 1.0, const double b = 1.0)
{
    vtkSmartPointer<vtkLineSource> lineSource = vtkSmartPointer<vtkLineSource>::New();
    lineSource->SetPoint1(point_a.x(), point_a.y(), point_a.z());
    lineSource->SetPoint2(point_b.x(), point_b.y(), point_b.z());

    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputConnection(lineSource->GetOutputPort());

    vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
    actor->SetMapper(mapper);
    actor->GetProperty()->SetColor(r, g, b);  // Red color for example
    actor->GetProperty()->SetLineWidth(1);

    return actor;
}

int main(int argc, char* argv[])
{
    for (int i = 0; i < argc; ++i) 
    {
        std::cout << argv[i] << std::endl;
    }

    vtkNew<vtkXMLPolyDataReader> vtp_reader;
    vtp_reader->SetFileName("arch.vtp");
    vtp_reader->Update();
    PolyData = vtp_reader->GetOutput();
    std::cout << "arch.vtp" << std::endl;
    vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    mapper->SetInputData(PolyData);
    mapper->Update();
    PolyDataActor = vtkSmartPointer<vtkActor>::New();
    PolyDataActor->SetMapper(mapper);
    PolyDataActor->GetProperty()->SetColor(1, 1, 1);
    PolyDataActor->GetProperty()->SetAmbient(0.5);
    PolyDataActor->GetProperty()->SetSpecularPower(100);
    PolyDataActor->GetProperty()->SetSpecular(0.5);
    PolyDataActor->GetProperty()->SetDiffuse(0.5);
    PolyDataActor->GetProperty()->EdgeVisibilityOff();
    PolyDataActor->PickableOn();
    Renderer->AddActor(PolyDataActor);

    vtkSmartPointer<vtkPLYReader> teeth_reader = vtkSmartPointer<vtkPLYReader>::New();
    teeth_reader->SetFileName("contacted.ply");
    teeth_reader->Update();
    TeethPolydata = teeth_reader->GetOutput();

    PolyDataToSurfaceMesh(TeethPolydata, teeth_sm);

    vtkSmartPointer<vtkPolyDataMapper> teeth_mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
    teeth_mapper->SetInputData(TeethPolydata);
    teeth_mapper->Update();
    TeethPolyDataActor = vtkSmartPointer<vtkActor>::New();
    TeethPolyDataActor->SetMapper(teeth_mapper);
    TeethPolyDataActor->GetProperty()->SetColor(0, 1, 1);
    TeethPolyDataActor->GetProperty()->SetAmbient(0.5);
    TeethPolyDataActor->GetProperty()->SetSpecularPower(100);
    TeethPolyDataActor->GetProperty()->SetSpecular(0.5);
    TeethPolyDataActor->GetProperty()->SetDiffuse(0.5);
    TeethPolyDataActor->GetProperty()->SetOpacity(0.8);
    //TeethPolyDataActor->GetProperty()->EdgeVisibilityOn();
    //TeethPolyDataActor->PickableOn();
    Renderer->AddActor(TeethPolyDataActor);

    Renderer->SetBackground(0.41, 0.41, 0.41);
    RenderWindow->AddRenderer(Renderer);
    RenderWindow->SetSize(1600, 900);
    RenderWindow->Render();
    Renderer->GetActiveCamera()->SetParallelProjection(1);
    vtkSmartPointer<DesignInteractorStyle> InteractorStyle = vtkSmartPointer<DesignInteractorStyle>::New();
    vtkSmartPointer<vtkRenderWindowInteractor> RenderWindowInteractor = vtkSmartPointer<vtkRenderWindowInteractor>::New();
    RenderWindowInteractor->SetInteractorStyle(InteractorStyle);
    RenderWindowInteractor->SetRenderWindow(RenderWindow);
    vtkNew<vtkCallbackCommand> LeftButtonPressCallback;
    LeftButtonPressCallback->SetCallback(LeftPress);
    LeftButtonPressCallback->InitializeObjectBase();
    vtkNew<vtkCallbackCommand> MouseMoveCallback;
    MouseMoveCallback->SetCallback(MouseMove);
    MouseMoveCallback->InitializeObjectBase();
    vtkNew<vtkCallbackCommand> LeftButtonReleaseCallback;
    LeftButtonReleaseCallback->SetCallback(LeftRelease);
    LeftButtonReleaseCallback->InitializeObjectBase();
    RenderWindowInteractor->AddObserver(vtkCommand::LeftButtonPressEvent, LeftButtonPressCallback);
    RenderWindowInteractor->AddObserver(vtkCommand::MouseMoveEvent, MouseMoveCallback);
    RenderWindowInteractor->AddObserver(vtkCommand::LeftButtonReleaseEvent, LeftButtonReleaseCallback);
    RenderWindowInteractor->Start();
}
