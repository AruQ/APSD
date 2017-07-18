cmake_minimum_required(VERSION 3.3)
project(Project_Sciddica)

find_package(MPI REQUIRED)
include_directories(${MPI_INCLUDE_PATH})

SET(CMAKE_C_COMPILER mpicc)
SET(CMAKE_CXX_COMPILER mpicxx)

set(SOURCE_FILES mainSciddica.cpp Reader.h)

add_executable(sciddica ${SOURCE_FILES})

TARGET_LINK_LIBRARIES(sciddica gomp)
