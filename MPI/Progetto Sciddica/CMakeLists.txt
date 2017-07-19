cmake_minimum_required(VERSION 2.8)
project(ColataLavica)

find_package(PkgConfig REQUIRED)
find_package(MPI REQUIRED)
pkg_search_module(GLFW REQUIRED glfw3)
pkg_search_module(ASSIMP REQUIRED assimp)
find_package(GLEW REQUIRED)
find_package(OpenGL  REQUIRED)

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++11")
include_directories(
  ${GLFW_INCLUDE_DIRS}
  ${GLEW_INCLUDE_DIRS}
  ${OPENGL_INCLUDE_DIRS}
  
include_directories(${MPI_INCLUDE_PATH})

  )

aux_source_directory(. SRC_LIST)

if(NOT MSVC)
  find_library(M_LIB m)
else() # not needed on Windows
  set(M_LIB "")
endif()

file(GLOB RES_FILE shaders/* *.frag *.vert *.vs *.glsl *.h)
add_executable(${PROJECT_NAME} ${SRC_LIST} ${RES_FILE})
target_link_libraries(
  ${PROJECT_NAME}
  ${GLFW_LIBRARIES}
  ${ASSIMP_LIBRARIES}
  ${GLEW_LIBRARIES}
  ${OPENGL_LIBRARIES}
  ${M_LIB}
  libSOIL.a
  )
TARGET_LINK_LIBRARIES(${PROJECT_NAME} gomp)
