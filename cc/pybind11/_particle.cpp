#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "particle.h"

namespace py = pybind11;

PYBIND11_MODULE(_particle, m){
    py::class_<Particle>(m, "Particle")
    .def(py::init<std::vector<float>,float>())
    .def(py::init<std::string,float>())
    
    .def_readonly("relPoints", &Particle::relPoints)
    .def_readonly("comPos", &Particle::comPos)
    .def_readonly("comVel", &Particle::comVel)
    .def_readonly("orient", &Particle::orient)
    .def_readonly("omega", &Particle::omega)
    
    .def("setComPos", &Particle::setComPos)
    .def("setComVel", &Particle::setComVel)
    .def("setOrient", &Particle::setOrient)
    .def("setOmega", &Particle::setOmega)
    
    .def("rotationalPotentialEnergy", &Particle::rotationalPotentialEnergy)
    .def("rotationalKineticEnergy", &Particle::rotationalKineticEnergy);
    //m.doc() = "";

}
