#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include "quat.h"

namespace py = pybind11;

PYBIND11_MODULE(_quat, m){
    auto m_quat = m.def_submodule("quat", "Submodule for quaternion operations");
    m_quat.def("multiply", py::overload_cast<Eigen::Vector4f,
                                             Eigen::Vector4f>(&multiply));
    m_quat.def("multiplyVecQuat", py::overload_cast<std::vector<float>,
                                                    std::vector<float>>(&multiplyVecQuat));
    m_quat.def("conj", py::overload_cast<std::vector<float>>(&conj));
    m_quat.def("vecRotate", py::overload_cast<std::vector<float>,
                                                    std::vector<float>>(&vecRotate));
    m.doc() = "";

}
