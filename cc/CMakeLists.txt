cmake_minimum_required(VERSION 3.24)
project("AggCol")

cmake_policy(SET CMP0135 NEW)
set(CMAKE_CXX_STANDARD 17)
enable_language(CXX)

set(CMAKE_INSTALL_PREFIX ${CMAKE_SOURCE_DIR}/../install)
set(CMAKE_INSTALL_RPATH ${CMAKE_INSTALL_PREFIX}/lib)
set(CMAKE_CXX_FLAGS "-03 -Wall")

include(FetchContent)

# get fmt library
FetchContent_Declare(
    fmt
    GIT_REPOSITORY https://github.com/fmtlib/fmt.git
    GIT_TAG v11.0.2
    )
    
# set external libs
set(EXTERNAL_LIBS fmt::fmt)

# includes
include_directories(${CMAKE_SOURCE_DIR}/include)

# project dirs
add_subdirectory(src)

if(PYTHON)
    add_subdirectory(pybind11)
endif()
