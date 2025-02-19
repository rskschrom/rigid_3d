#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "particle.h"

namespace py = pybind11;

PYBIND11_MODULE(_particle, m){
    py::class_<Particle>(m, "Particle")
    .def(py::init<std::vector<float>,float>())
    .def(py::init<std::string,float>())
    
    .def_readonly("points", &Particle::points)
    .def_readonly("relPoints", &Particle::relPoints)
        
    .def("inertiaMomentTensor", &Particle::inertiaMomentTensor)
    .def("getMatInerm", &Particle::getMatInerm)
    .def("setMatInerm", &Particle::setMatInerm);
    //m.doc() = "";

}
