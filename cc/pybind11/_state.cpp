#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "state.h"

namespace py = pybind11;

PYBIND11_MODULE(_state, m){
    py::class_<State>(m, "State")
    .def(py::init<Eigen::Matrix3f,std::vector<float>,std::vector<float>,
                  std::vector<float>,std::vector<float>>())
    .def(py::init<Eigen::Matrix3f>())
    
    .def_readonly("comPos", &State::comPos)
    .def_readonly("comVel", &State::comVel)
    .def_readonly("orient", &State::orient)
    .def_readonly("omega", &State::omega)
    
    .def("setComPos", &State::setComPos)
    .def("setComVel", &State::setComVel)
    .def("setOrient", &State::setOrient)
    .def("setOmega", &State::setOmega)
    
    .def("rotationalPotentialEnergy", &State::rotationalPotentialEnergy)
    .def("rotationalKineticEnergy", &State::rotationalKineticEnergy);
    //m.doc() = "";

}
