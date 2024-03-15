#if 1
#pragma once
#include <iostream>
#include <algorithm>
#include <cmath>
#include <vtkOBBTree.h>
#include <vtkLineSource.h>
#include <vtkTubeFilter.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkFeatureEdges.h>
#include <vtkSTLReader.h>
#include <vtkPLYReader.h>
#include <vtkOBJReader.h>
#include <vtkSTLWriter.h>
#include <vtkPLYWriter.h>
#include <vtkNamedColors.h>
#include <vtkSmartPointer.h>
#include <vtkNew.h>
#include <vtkPolyData.h>
#include <vtkPolyDataMapper.h>
#include <vtkActor.h>
#include <vtkCamera.h>
#include <vtkProperty.h>
#include <vtkRenderer.h>
#include <vtkRenderWindow.h>
#include <vtkRenderWindowInteractor.h>
#include <vtkInteractorStyleTrackballCamera.h>
#include <vtkPointData.h>
#include <vtkUnsignedCharArray.h>
#include <vtkMath.h>
#include <vtkArrowSource.h>
#include <vtkTransformPolyDataFilter.h>
#include<vtkCallbackCommand.h>
#include<vtkCellPicker.h>
#include<vtkPolyDataNormals.h>
#include<vtkSphereSource.h>
#include <vtkPolyLine.h>
#include<vtkParametricSpline.h>
#include<vtkParametricFunctionSource.h>
#include<vtkCellArray.h>
#include<vtkTriangle.h>
#include<vtkMath.h>
#include<vtkCellData.h>
#include<vector>
#include <vtkAutoInit.h>

#include <CGAL/Polygon_mesh_processing/triangulate_hole.h>
#include <CGAL/Polygon_mesh_processing/corefinement.h>
#include <CGAL/IO/write_ply_points.h>
VTK_MODULE_INIT(vtkRenderingOpenGL2)
VTK_MODULE_INIT(vtkInteractionStyle)
VTK_MODULE_INIT(vtkRenderingFreeType)
#include <CGAL/Simple_cartesian.h>
#include <CGAL/Exact_predicates_exact_constructions_kernel.h>
#include <CGAL/Surface_mesh.h>
#include <CGAL/Surface_mesh_parameterization/IO/File_off.h>
#include <CGAL/Surface_mesh_parameterization/parameterize.h>
#include <CGAL/Polygon_mesh_processing/measure.h>
#include <CGAL/Surface_mesh_deformation.h>
#include <CGAL/Polygon_mesh_processing/remesh.h>
#include <CGAL/Segment_2.h>
#include <CGAL/intersections.h>
#include <CGAL/Polygon_mesh_processing/compute_normal.h>
#include <CGAL/Polygon_mesh_processing/extrude.h>
#include <CGAL/boost/graph/Euler_operations.h>
#include <iostream>
#include <fstream>
#include <string>
#include <omp.h>
typedef CGAL::Simple_cartesian<double>       Kernel;
typedef Kernel::Point_2                         Point_2;
typedef Kernel::Point_3                         Point_3;
typedef Kernel::Vector_3						Vector_3;
typedef CGAL::Surface_mesh<Kernel::Point_3>     SurfaceMesh;
typedef boost::graph_traits<SurfaceMesh>::vertex_descriptor     vertex_descriptor;
typedef boost::graph_traits<SurfaceMesh>::halfedge_descriptor   halfedge_descriptor;
typedef boost::graph_traits<SurfaceMesh>::vertex_iterator        vertex_iterator;
typedef boost::graph_traits<SurfaceMesh>::face_descriptor       face_descriptor;
typedef boost::graph_traits<SurfaceMesh>::edge_descriptor       edge_descriptor;
typedef CGAL::Surface_mesh_deformation<SurfaceMesh, CGAL::Default, CGAL::Default, CGAL::SRE_ARAP> Surface_mesh_deformation;
typedef boost::property_map<SurfaceMesh, CGAL::vertex_point_t>::type  VPMap;
typedef SurfaceMesh::template Property_map<vertex_descriptor, Vector_3> VNMap;
typedef SurfaceMesh::template Property_map<vertex_descriptor, double> VLMap;

namespace SMP = CGAL::Surface_mesh_parameterization;
namespace PMP = CGAL::Polygon_mesh_processing;
SurfaceMesh sm;
SurfaceMesh cut_sm;
SurfaceMesh teeth_sm;
SurfaceMesh bite_sm;
std::map<unsigned int, face_descriptor>FaceMap;
std::map<unsigned int, vertex_descriptor>VertexMap;
std::map<unsigned int, halfedge_descriptor>HEMap;

std::vector<int> vVisit;
std::vector<int> eVisit;
int eIndex = 0;
int vIndex = 0;
SurfaceMesh::Property_map<vertex_descriptor, Point_2> uv_map = sm.add_property_map<vertex_descriptor, Point_2>("h:uv").first;


typedef Kernel::Vector_3                                     Vector_3;
typedef boost::property_map<SurfaceMesh, CGAL::vertex_point_t>::type  VPMap;
typedef SurfaceMesh::template Property_map<vertex_descriptor, Vector_3> VNMap;
typedef SurfaceMesh::template Property_map<vertex_descriptor, double> VLMap;

struct Bottom
{
	Bottom(VPMap pmap, VNMap nmap, VLMap lmap)
		: pmap(pmap), nmap(nmap), lmap(lmap)
	{}
	void operator()(const vertex_descriptor& vin, const vertex_descriptor vout) const
	{
		put(pmap, vout, get(pmap, vout) - get(lmap, vin) * get(nmap, vin));
	}
	VPMap pmap;
	VNMap nmap;
	VLMap lmap;

};
struct Top
{
	Top(VPMap pmap, VNMap nmap, VLMap lmap)
		: pmap(pmap), nmap(nmap), lmap(lmap)
	{}
	void operator()(const vertex_descriptor& vin, const vertex_descriptor vout) const
	{
		put(pmap, vout, get(pmap, vout) + get(lmap, vin) * get(nmap, vin));
	}
	VPMap pmap;
	VNMap nmap;
	VLMap lmap;

};

vtkSmartPointer<vtkPolyData> CutPolydata = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> TeethPolydata = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyData> BitePolyData = vtkSmartPointer<vtkPolyData>::New();

vtkSmartPointer<vtkPolyData> PolyData = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkRenderWindow> RenderWindow = vtkSmartPointer<vtkRenderWindow>::New();
vtkSmartPointer<vtkRenderer> Renderer = vtkSmartPointer<vtkRenderer>::New();
vtkSmartPointer<vtkOBBTree> PolyData_OBBTree = vtkSmartPointer<vtkOBBTree>::New();
vtkSmartPointer<vtkCellPicker> CellPicker = vtkSmartPointer<vtkCellPicker>::New();
vtkSmartPointer<vtkPolyData> Spline = vtkSmartPointer<vtkPolyData>::New();
vtkSmartPointer<vtkPolyLine> PolyLine = vtkSmartPointer<vtkPolyLine>::New();
vtkSmartPointer<vtkPoints> Points = vtkSmartPointer<vtkPoints>::New();
vtkSmartPointer<vtkActor> SplineActor = vtkSmartPointer<vtkActor>::New();
SurfaceMesh teeth;

struct MeshPoint
{
	unsigned int TriId;
	Point_2 uv;
	double xyz[3];
};

std::vector<MeshPoint> MeshSpline;
class DesignInteractorStyle : public vtkInteractorStyleTrackballCamera
{
public:
	static DesignInteractorStyle* New()
	{
		return new DesignInteractorStyle;
	}
	vtkTypeMacro(DesignInteractorStyle, vtkInteractorStyleTrackballCamera);

	DesignInteractorStyle() {}
	virtual ~DesignInteractorStyle() {}
	virtual void OnLeftButtonDown() {}
	virtual void OnLeftButtonUp() {}
	virtual void OnRightButtonDown() { this->StartRotate(); } // 避免vtk的GrabFocus接口占用交互命令
	virtual void OnRightButtonUp() { this->vtkInteractorStyleTrackballCamera::OnLeftButtonUp(); }
	virtual void OnMouseMove() { this->vtkInteractorStyleTrackballCamera::OnMouseMove(); }
	virtual void OnMouseWheelForward() { if (!this->Interactor->GetControlKey() && !this->Interactor->GetShiftKey()) this->vtkInteractorStyleTrackballCamera::OnMouseWheelForward(); }
	virtual void OnMouseWheelBackward() { if (!this->Interactor->GetControlKey() && !this->Interactor->GetShiftKey()) this->vtkInteractorStyleTrackballCamera::OnMouseWheelBackward(); }
};

vtkSmartPointer<vtkActor> MakeActor(vtkSmartPointer<vtkPolyData> polydata, double R, double G, double B)
{
	vtkSmartPointer<vtkActor> actor = vtkSmartPointer<vtkActor>::New();
	if (!polydata) return actor;
	vtkSmartPointer<vtkPolyDataNormals> vtkNormal = vtkSmartPointer<vtkPolyDataNormals>::New();
	vtkNormal->SetInputData(polydata);
	vtkNormal->SetComputePointNormals(1);
	vtkNormal->SetComputeCellNormals(1);
	vtkNormal->SetAutoOrientNormals(1);
	vtkNormal->SetSplitting(0);
	vtkNormal->FlipNormalsOn();
	vtkNormal->Update();
	vtkSmartPointer<vtkPolyDataMapper> mapper = vtkSmartPointer<vtkPolyDataMapper>::New();
	mapper->SetInputData(vtkNormal->GetOutput());
	mapper->Update();
	actor->SetMapper(mapper);
	actor->GetProperty()->SetColor(R, G, B);
	actor->GetProperty()->SetAmbient(0.5);
	actor->GetProperty()->SetSpecularPower(100);
	actor->GetProperty()->SetSpecular(0.5);
	actor->GetProperty()->SetDiffuse(0.5);
	actor->PickableOff();
	
	return actor;
}

void GetRayIntersection(vtkObject* caller, int& nTriId, double* pWorld, vtkSmartPointer<vtkActor> &actor)
{
	vtkSmartPointer<vtkRenderWindowInteractor>vtkInter = vtkSmartPointer<vtkRenderWindowInteractor>::New();
	vtkInter = vtkRenderWindowInteractor::SafeDownCast(caller);
	if (!vtkInter)
	{
		return;
	}
	int* pEvtPos = vtkInter->GetEventPosition();
	vtkInter->FindPokedRenderer(pEvtPos[0], pEvtPos[1]);
	vtkInter->SetPicker(CellPicker);
	vtkInter->GetPicker()->Pick(pEvtPos[0], pEvtPos[1], 0, Renderer.GetPointer());

	nTriId = CellPicker->GetCellId();
	CellPicker->GetPickPosition(pWorld);
	actor = CellPicker->GetActor();
}


bool GetBary(double v1[3], double v2[3], double v3[3], double vp[3], double& Bary1, double& Bary2, double& Bary3)
{
	double v02[3];
	v02[0] = v1[0] - v3[0];		v02[1] = v1[1] - v3[1];		v02[2] = v1[2] - v3[2];
	double v12[3];
	v12[0] = v2[0] - v3[0];		v12[1] = v2[1] - v3[1];		v12[2] = v2[2] - v3[2];
	double vp2[3];
	vp2[0] = vp[0] - v3[0];		vp2[1] = vp[1] - v3[1];		vp2[2] = vp[2] - v3[2];
	double d0202 = v02[0] * v02[0] + v02[1] * v02[1] + v02[2] * v02[2];
	double d0212 = v02[0] * v12[0] + v02[1] * v12[1] + v02[2] * v12[2];
	double d1212 = v12[0] * v12[0] + v12[1] * v12[1] + v12[2] * v12[2];
	double d02p2 = v02[0] * vp2[0] + v02[1] * vp2[1] + v02[2] * vp2[2];
	double d12p2 = v12[0] * vp2[0] + v12[1] * vp2[1] + v12[2] * vp2[2];
	double Det = d0202 * d1212 - d0212 * d0212;

	double dInvDet = static_cast<double>(1.0) / Det;

	Bary1 = (d1212 * d02p2 - d0212 * d12p2) / Det;
	Bary2 = (d0202 * d12p2 - d0212 * d02p2) / Det;
	Bary3 = static_cast<double>(1.0) - Bary1 - Bary2;

	return ((Bary1 >= 0) && (Bary2 >= 0) && (Bary1 + Bary2 < static_cast<double>(1.0)));
}

Point_2 EdgeCross(Point_2 p1, Point_2 p2, Point_2 p3, Point_2 p4)
{
	Kernel::Segment_2 seg1(p1, p2);
	Kernel::Segment_2 seg2(p3, p4);
	auto ans = CGAL::intersection(seg1, seg2);
	if (ans)
	{
		if (const Kernel::Point_2* s = boost::get<Kernel::Point_2>(&*ans))
		{
			Point_2 cross(s->x(), s->y());
			return cross;
		}
		//std::cout << *ans <<" 0" << std::endl;

	}
	else
	{
		return Point_2(0, 0);
	}
}

bool inTriangle(int TriId, Point_2 vp)
{
	unsigned int vid0 = PolyData->GetCell(TriId)->GetPointIds()->GetId(0);
	unsigned int vid1 = PolyData->GetCell(TriId)->GetPointIds()->GetId(1);
	unsigned int vid2 = PolyData->GetCell(TriId)->GetPointIds()->GetId(2);
	auto v1 = uv_map[VertexMap[vid0]];
	auto v2 = uv_map[VertexMap[vid1]];
	auto v3 = uv_map[VertexMap[vid2]];
	double v02[3];
	v02[0] = v1[0] - v3[0];		v02[1] = v1[1] - v3[1];
	double v12[3];
	v12[0] = v2[0] - v3[0];		v12[1] = v2[1] - v3[1];
	double vp2[3];
	vp2[0] = vp[0] - v3[0];		vp2[1] = vp[1] - v3[1];
	double d0202 = v02[0] * v02[0] + v02[1] * v02[1];
	double d0212 = v02[0] * v12[0] + v02[1] * v12[1];
	double d1212 = v12[0] * v12[0] + v12[1] * v12[1];
	double d02p2 = v02[0] * vp2[0] + v02[1] * vp2[1];
	double d12p2 = v12[0] * vp2[0] + v12[1] * vp2[1];
	double Det = d0202 * d1212 - d0212 * d0212;

	double dInvDet = static_cast<double>(1.0) / Det;

	double Bary1 = (d1212 * d02p2 - d0212 * d12p2) / Det;
	double Bary2 = (d0202 * d12p2 - d0212 * d02p2) / Det;
	double Bary3 = static_cast<double>(1.0) - Bary1 - Bary2;

	return ((Bary1 >= 0) && (Bary2 >= 0) && (Bary1 + Bary2 <= static_cast<double>(1.0)));
}

void UV2XYZ(int TriId, Point_2 uv,double ans[3])
{
	unsigned int vid0 = PolyData->GetCell(TriId)->GetPointIds()->GetId(0);
	unsigned int vid1 = PolyData->GetCell(TriId)->GetPointIds()->GetId(1);
	unsigned int vid2 = PolyData->GetCell(TriId)->GetPointIds()->GetId(2);
	auto v1 = uv_map[VertexMap[vid0]];
	auto v2 = uv_map[VertexMap[vid1]];
	auto v3 = uv_map[VertexMap[vid2]];
	double v02[3];
	v02[0] = v1[0] - v3[0];		v02[1] = v1[1] - v3[1];
	double v12[3];
	v12[0] = v2[0] - v3[0];		v12[1] = v2[1] - v3[1];
	double vp2[3];
	vp2[0] = uv[0] - v3[0];		vp2[1] = uv[1] - v3[1];
	double d0202 = v02[0] * v02[0] + v02[1] * v02[1];
	double d0212 = v02[0] * v12[0] + v02[1] * v12[1];
	double d1212 = v12[0] * v12[0] + v12[1] * v12[1];
	double d02p2 = v02[0] * vp2[0] + v02[1] * vp2[1];
	double d12p2 = v12[0] * vp2[0] + v12[1] * vp2[1];
	double Det = d0202 * d1212 - d0212 * d0212;

	double dInvDet = static_cast<double>(1.0) / Det;

	double Bary1 = (d1212 * d02p2 - d0212 * d12p2) / Det;
	double Bary2 = (d0202 * d12p2 - d0212 * d02p2) / Det;
	double Bary3 = static_cast<double>(1.0) - Bary1 - Bary2;
	double v3d0[3], v3d1[3], v3d2[3];
	PolyData->GetCell(TriId)->GetPoints()->GetPoint(0, v3d0);
	PolyData->GetCell(TriId)->GetPoints()->GetPoint(1, v3d1);
	PolyData->GetCell(TriId)->GetPoints()->GetPoint(2, v3d2);
	ans[0] = v3d0[0] * Bary1 + v3d1[0] * Bary2 + v3d2[0] * Bary3;
	ans[1] = v3d0[1] * Bary1 + v3d1[1] * Bary2 + v3d2[1] * Bary3;
	ans[2] = v3d0[2] * Bary1 + v3d1[2] * Bary2 + v3d2[2] * Bary3;
}

halfedge_descriptor FindAllCrosses(halfedge_descriptor he, Point_2 uv1, Point_2 uv2, bool first)
{
	if (inTriangle(sm.face(he).idx(), uv2))
	{
		MeshPoint p;
		double pt[3];
		p.TriId = sm.face(he).idx();
		p.uv = uv1;
		UV2XYZ(p.TriId, p.uv, pt);
		p.xyz[0] = pt[0]; p.xyz[1] = pt[1]; p.xyz[2] = pt[2];
		MeshSpline.push_back(p);

		MeshPoint p1;
		double pt1[3];
		p1.TriId = sm.face(he).idx();
		p1.uv = uv2;
		UV2XYZ(p1.TriId, p1.uv, pt1);
		p1.xyz[0] = pt1[0]; p1.xyz[1] = pt1[1]; p1.xyz[2] = pt1[2];
		MeshSpline.push_back(p1);
		//cout << "end" << endl;
		return he;
	}
	auto he0 = he;
	Point_2 crossVer;
	if (first)
	{
		crossVer = EdgeCross(uv1, uv2, uv_map[sm.source(he0)], uv_map[sm.target(he0)]);
		if (crossVer != Point_2(0, 0))
		{
			MeshPoint p;
			double pt[3];
			p.TriId = sm.face(he).idx();
			p.uv = uv1;
			UV2XYZ(p.TriId, p.uv, pt);
			p.xyz[0] = pt[0]; p.xyz[1] = pt[1]; p.xyz[2] = pt[2];
			MeshSpline.push_back(p);

			return FindAllCrosses(sm.opposite(he0), crossVer, uv2, 0);
		}
	}
	auto he1 = sm.next(he0);
	crossVer = EdgeCross(uv1, uv2, uv_map[sm.source(he1)], uv_map[sm.target(he1)]);
	if (crossVer != Point_2(0, 0))
	{
		MeshPoint p;
		double pt[3];
		p.TriId = sm.face(he).idx();
		p.uv = uv1;
		UV2XYZ(p.TriId, p.uv, pt);
		p.xyz[0] = pt[0]; p.xyz[1] = pt[1]; p.xyz[2] = pt[2];
		MeshSpline.push_back(p);
		//cout << crossVer << endl;
		
		return FindAllCrosses(sm.opposite(he1), crossVer, uv2, 0);
	}
	auto he2 = sm.next(he1);
	crossVer = EdgeCross(uv1, uv2, uv_map[sm.source(he2)], uv_map[sm.target(he2)]);
	if (crossVer != Point_2(0, 0))
	{
		MeshPoint p;
		double pt[3];
		p.TriId = sm.face(he).idx();
		p.uv = uv1;
		UV2XYZ(p.TriId, p.uv, pt);
		p.xyz[0] = pt[0]; p.xyz[1] = pt[1]; p.xyz[2] = pt[2];
		MeshSpline.push_back(p);
		//cout << crossVer << endl;
		return FindAllCrosses(sm.opposite(he2), crossVer, uv2, 0);
	}
	//cout << uv1 << endl;
	//cout << uv2 << endl;
	//cout << "error" << endl;
	//return he;
}

//CGAL转换vtkPolydata
vtkSmartPointer<vtkPolyData> CGAL_Surface_Mesh2VTK_PolyData(const SurfaceMesh& pmesh)
{
	vtkSmartPointer<vtkPoints>  Points = vtkSmartPointer<vtkPoints>::New();
	Points->SetNumberOfPoints(pmesh.number_of_vertices());
	//vtkSmartPointer<vtkCellArray> TrianglePolys = vtkSmartPointer<vtkCellArray>::New();
	//TrianglePolys->SetNumberOfCells(pmesh.number_of_faces());
	vtkSmartPointer<vtkIdTypeArray> connectivity = vtkSmartPointer<vtkIdTypeArray>::New();
	connectivity->SetNumberOfComponents(4);
	connectivity->SetNumberOfTuples(pmesh.number_of_faces());
	vtkSmartPointer<vtkPolyData> Polydata = vtkSmartPointer<vtkPolyData>::New();
	//记录点描述符和index的映射表VertexDescriptorIndexMap
	std::map<vertex_descriptor, int> VertexDescriptorIndexMap;
	int VertexIndex = 0;
	int TriIndex = 0;

	//插入点(遍历所有有效点)

	for (vertex_descriptor& v : CGAL::vertices(pmesh))
	{
		//CGAL::get(CGAL::vertex_point, pmesh)返回一张<顶点描述符,点坐标>的映射表  传入描述符v得到点坐标p
		const boost::property_map_value<SurfaceMesh, CGAL::vertex_point_t>::type& p = get(CGAL::get(CGAL::vertex_point, pmesh), v);
		//Points->InsertNextPoint(CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
		Points->InsertPoint(VertexIndex,CGAL::to_double(p.x()), CGAL::to_double(p.y()), CGAL::to_double(p.z()));
		//map中插入<描述符,索引>
		VertexDescriptorIndexMap[v] = VertexIndex++;
	}

	//插入面
	for (face_descriptor& f : CGAL::faces(pmesh))
	{
		//TriIndex = 0;
		//vtkSmartPointer<vtkTriangle> triangle = vtkSmartPointer<vtkTriangle>::New();
		//找到pmesh中指向面描述符f的半边 再调用halfedges_around_face(f->next->next逆时针遍历face)
		//CGAL::target(h, pmesh):半边描述符h指向的点描述符
		//for (halfedge_descriptor h : CGAL::halfedges_around_face(CGAL::halfedge(f, pmesh), pmesh))
		//	triangle->GetPointIds()->SetId(TriIndex++, VertexDescriptorIndexMap[CGAL::target(h, pmesh)]);
		int n = 0;
		int ids[3];
		for (halfedge_descriptor &h : CGAL::halfedges_around_face(CGAL::halfedge(f, pmesh), pmesh))
		{
			ids[n++] = VertexDescriptorIndexMap[CGAL::target(h, pmesh)];
		}
		connectivity->SetTuple4(TriIndex++, 3, ids[0], ids[1], ids[2]);
		//TrianglePolys->InsertNextCell(triangle);
	}
	vtkSmartPointer<vtkCellArray> cells = vtkSmartPointer<vtkCellArray>::New();
	cells->SetCells(pmesh.number_of_faces(), connectivity);
	Polydata->SetPoints(Points);
	Polydata->SetPolys(cells);
	return Polydata;
}

struct Facet_with_normal_pmap
{
	template <class Facet>
	struct Bind
	{
		typedef typename Facet::Normal value_type;
		typedef typename Facet::ID ID;
		typedef boost::read_write_property_map_tag category;
		typedef const value_type& reference;
	};

	template <class Facet_with_id>
	Facet_with_id get(const typename Facet_with_id::ID id) const
	{
		return Facet_with_id(id, Facet_with_id::facet(id).normal());
	}

	template <class Facet_with_id>
	void put(const typename Facet_with_id::ID id, const typename Facet_with_id::value_type& value) const
	{
		Facet_with_id::facet(id).set_normal(value);
	}
};
#endif