#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <pybind11/eigen.h>
#include "mesh_particle.h"

namespace py = pybind11;

PYBIND11_MODULE(_mesh_particle, m){
    py::class_<MeshParticle>(m, "MeshParticle")
    .def(py::init<Eigen::MatrixX3f,Eigen::MatrixX3i,float>())
    
    .def_readonly("vertices", &MeshParticle::vertices)
    .def_readonly("faces", &MeshParticle::faces)
    .def_readonly("faceAreas", &MeshParticle::faceAreas)    
    .def_readonly("faceNorms", &MeshParticle::faceNorms)
        
    .def("totalMass", &MeshParticle::totalMass)
    .def("triAlp2BetIntegral", &MeshParticle::triAlp2BetIntegral)
    .def("calculateFaceAreasNorms", &MeshParticle::calculateFaceAreasNorms)
    //.def("inertiaMomentTensor", &MeshParticle::inertiaMomentTensor)
    .def("getMatInerm", &MeshParticle::getMatInerm)
    .def("setMatInerm", &MeshParticle::setMatInerm);
    //m.doc() = "";

}
