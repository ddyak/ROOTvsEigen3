cmake_minimum_required (VERSION 3.1)
project (benchmarks)

#//g++ main.cc -std=c++11 -lbenchmark -lpthread -O2 -o benchmark

set(CMAKE_CXX_STANDARD 14)
set(CMAKE_CXX_STANDARD_REQUIRED ON)
set(CMAKE_CXX_EXTENSIONS OFF)
set(CMAKE_BUILD_TYPE Release)

find_package(Threads)
find_package(benchmark REQUIRED)

find_package(ROOT REQUIRED COMPONENTS EG)
find_package(Eigen3 REQUIRED)

link_directories(
   /home/ddyak/sctanalysis/utility/build/
)

include_directories(
    include
    /home/ddyak/sctanalysis/utility/include
    ${EIGEN3_INCLUDE_DIR}
)

file(GLOB SOURCES src/*.cc)

add_executable(${PROJECT_NAME} src/main.cc)

# add_library(${PROJECT_NAME} SHARED ${SOURCES})
target_link_libraries(${PROJECT_NAME} ${ROOT_LIBRARIES} libbenchmark.so ${CMAKE_THREAD_LIBS_INIT} libanaUtil.so)