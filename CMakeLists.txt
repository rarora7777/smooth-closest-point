cmake_minimum_required(VERSION 3.13)
project(PhongProjection)

set(CMAKE_MODULE_PATH ${CMAKE_MODULE_PATH} ${CMAKE_CURRENT_SOURCE_DIR}/cmake)

find_package(OpenMesh REQUIRED)
find_package(Eigen3 3.3 REQUIRED NO_MODULE)
find_package(LIBIGL REQUIRED)

# Add your project files

add_executable(decimator decimator/decimator.cc)
target_link_directories(decimator PUBLIC ${OPENMESH_LIBRARIES_DIR})
target_link_libraries(decimator ${OPENMESH_LIBRARIES})
target_include_directories(decimator PUBLIC ${OPENMESH_INCLUDE_DIRS})
set_target_properties(decimator PROPERTIES RUNTIME_OUTPUT_DIRECTORY_RELEASE ${CMAKE_CURRENT_SOURCE_DIR}/decimator)
set_target_properties(decimator PROPERTIES RUNTIME_OUTPUT_DIRECTORY_DEBUG ${CMAKE_CURRENT_SOURCE_DIR}/decimator)
set_target_properties(decimator PROPERTIES RUNTIME_OUTPUT_DIRECTORY ${CMAKE_CURRENT_SOURCE_DIR}/decimator)


add_library(Phong SHARED src/phong/Phong.cpp src/phong/TrianglePhong.cpp)
target_include_directories(Phong PUBLIC ${CMAKE_CURRENT_SOURCE_DIR}/src ${CMAKE_CURRENT_SOURCE_DIR}/src/utils ${LIBIGL_INCLUDE_DIR})
target_link_libraries(Phong Eigen3::Eigen igl::core)