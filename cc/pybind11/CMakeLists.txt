include_directories(${CMAKE_SOURCE_DIR/src})
FetchContent_Declare(
    pybind11
    GIT_REPOSITORY https://github.com/pybind/pybind11
    GIT_TAG v2.12
    EXCLUDE_FROM_ALL
)
FetchContent_MakeAvailable(pybind11)

find_package(Python COMPONENTS Interpreter Development)

set(CMAKE_INSTALL_RPATH_USE_LINK_PATH)

# loop over modules to build
FILE(GLOB modules _*.cpp)

foreach(mod ${modules})
    cmake_path(GET mod FILENAME mod_name)
    string(REPLACE ".cpp" "" mod_name ${mod_name})
    pybind11_add_module(${mod_name} MODULE ${mod})
    target_link_libraries(${mod_name} PRIVATE aggcol pybind11::pybind11)
    install(TARGETS ${mod_name} LIBRARY DESTINATION ${CMAKE_SOURCE_DIR}/../py/src/pybind11_lib)
endforeach()
