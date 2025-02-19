#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "freefallsimulation.h"

namespace py = pybind11;

PYBIND11_MODULE(_freefallsimulation, m){
    py::class_<FreeFallSimulation>(m, "FreeFallSimulation")
    .def(py::init<State,float,float,float>())
    .def(py::init<State,int,float,float,float,float,float>())
    
    .def_readonly("posHistory", &FreeFallSimulation::posHistory)
    .def_readonly("orientHistory", &FreeFallSimulation::orientHistory)
    .def_readonly("st", &FreeFallSimulation::st)
        
    .def("getSimStep", &FreeFallSimulation::getSimStep)
    
    .def("evolveMotionInertial", &FreeFallSimulation::evolveMotionInertial)
    .def("evolveMotionBuoyancy", &FreeFallSimulation::evolveMotionBuoyancy);
    //m.doc() = "";

}
