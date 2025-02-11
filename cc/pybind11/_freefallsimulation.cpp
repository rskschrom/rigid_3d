#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "freefallsimulation.h"
#include "particle.h"

namespace py = pybind11;

PYBIND11_MODULE(_freefallsimulation, m){
    py::class_<FreeFallSimulation>(m, "FreeFallSimulation")
    .def(py::init<Particle,float,float,float>())
    .def(py::init<Particle,int,float,float,float,float,float>())
    
    .def_readonly("posHistory", &FreeFallSimulation::posHistory)
    .def_readonly("orientHistory", &FreeFallSimulation::orientHistory)
    .def_readonly("par", &FreeFallSimulation::par)
        
    .def("getSimStep", &FreeFallSimulation::getSimStep)
    
    .def("evolveMotionInertial", &FreeFallSimulation::evolveMotionInertial)
    .def("evolveMotionBuoyancy", &FreeFallSimulation::evolveMotionBuoyancy);
    //m.doc() = "";

}
