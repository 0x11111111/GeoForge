# Set the minimum required version of CMake
cmake_minimum_required(VERSION 3.10)

# Set the project name and version
project(GeoForge VERSION 1.0)

# Set the C++ standard
set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED True)

# Find the Qt library
find_package(Qt5 COMPONENTS Widgets Core Gui OpenGL REQUIRED)

# Find the VTK library
find_package(VTK REQUIRED)
include(${VTK_USE_FILE})

# Find the CGAL library
find_package(CGAL REQUIRED COMPONENTS Core)
include(${CGAL_USE_FILE})

# Find the Eigen3 library
find_package(Eigen3 3.3 REQUIRED NO_MODULE)

# Set the optimization level to O2 for Release mode
if(MSVC)
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} /O2")
else()
  set(CMAKE_CXX_FLAGS_RELEASE "${CMAKE_CXX_FLAGS_RELEASE} -O2")
endif()

# Find all sources
file(GLOB_RECURSE SOURCE_FILES
    "${CMAKE_SOURCE_DIR}/src/*.cpp"
    "${CMAKE_SOURCE_DIR}/src/*.hpp"
    "${CMAKE_SOURCE_DIR}/src/*.h"
    "${CMAKE_SOURCE_DIR}/src/*.ui"
)

# Add an executable for your project
add_executable(${PROJECT_NAME} ${SOURCE_FILES})

# Link libraries to your project
target_link_libraries(${PROJECT_NAME}
  Qt5::Widgets
  Qt5::Core
  Qt5::Gui
  Qt5::OpenGL
  ${VTK_LIBRARIES}
  CGAL::CGAL CGAL::CGAL_Core
  Eigen3::Eigen
)

# If Libigl is enabled, link the Libigl libraries
if(WITH_LIBIGL)
  target_link_libraries(${PROJECT_NAME} igl::core igl::opengl igl::opengl_glfw)
endif()

# Set to generate debug information
if(MSVC)
  set(CMAKE_SHARED_LINKER_FLAGS_RELEASE "${CMAKE_SHARED_LINKER_FLAGS_RELEASE} /DEBUG")
  set(CMAKE_EXE_LINKER_FLAGS_RELEASE "${CMAKE_EXE_LINKER_FLAGS_RELEASE} /DEBUG")
endif()

if(MSVC)
  # Set stack size
  set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} /STACK:204857600")
  # Set subsystem
  set_target_properties(${PROJECT_NAME} PROPERTIES
      LINK_FLAGS_RELEASE "/SUBSYSTEM:CONSOLE"
  )

  add_compile_options(/Zi)
endif()

# Configure installation rules
install(TARGETS ${PROJECT_NAME} DESTINATION bin)