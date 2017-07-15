cmake_minimum_required(VERSION 2.8)

project(Project)

add_executable(mytest mainSciddica.cpp)

# Require MPI for this project:
find_package(MPI REQUIRED)
set(CMAKE_CXX_COMPILE_FLAGS ${CMAKE_CXX_COMPILE_FLAGS} ${MPI_COMPILE_FLAGS})
set(CMAKE_CXX_LINK_FLAGS ${CMAKE_CXX_LINK_FLAGS} ${MPI_LINK_FLAGS})
include_directories(MPI_INCLUDE_PATH)
target_link_libraries(mytest ${MPI_LIBRARIES}) 


# Add a test:
enable_testing()
add_test(SimpleTest ${MPIEXEC} ${MPIEXEC_NUMPROC_FLAG} 4 ${CMAKE_CURRENT_BINARY_DIR}/mytest)
