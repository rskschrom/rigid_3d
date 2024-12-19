cmake_minimum_required(VERSION 3.28)
project("AggCol")

cmake_policy(SET CMP0135 NEW)
set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../install)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_CXX_FLAGS "-O3 -Wall -fPIC")

include(FetchContent)

# get fmt library
#FetchContent_Declare(
#    fmt
#    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
#    GIT_TAG 11.0.2
#    )
    
# get the eigen library
FetchContent_Declare(
    Eigen
    GIT_REPOSITORY https://gitlab.com/libeigen/eigen.git
    GIT_TAG master
    EXCLUDE_FROM_ALL
    )

FetchContent_MakeAvailable(Eigen)  
#FetchContent_MakeAvailable(fmt)

# set external libs
#set(EXTERNAL_LIBS fmt::fmt Eigen3::Eigen)
set(EXTERNAL_LIBS Eigen3::Eigen)

# includes
include_directories(${CMAKE_SOURCE_DIR}/include)

# project dirs
add_subdirectory(src)
add_subdirectory(examples)

if(PYTHON)
    add_subdirectory(pybind11)
endif()
