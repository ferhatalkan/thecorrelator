cmake_minimum_required(VERSION 3.10)
project(thecorrelator VERSION 1.0 LANGUAGES CXX)

set(CMAKE_CXX_STANDARD 11)
set(CMAKE_CXX_STANDARD_REQUIRED ON)

# OpenMP
find_package(OpenMP REQUIRED)

# GSL
find_package(GSL REQUIRED)

# zlib
find_package(ZLIB REQUIRED)

add_executable(correlate correlator_parallel.cpp)

target_link_libraries(correlate
    PRIVATE
    OpenMP::OpenMP_CXX
    GSL::gsl
    GSL::gslcblas
    ZLIB::ZLIB
)
