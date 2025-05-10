#include <polyfem/mesh/SlimSmooth.hpp>
#include <polyfem/solver/AdjointTools.hpp>
#include <polyfem/utils/MatrixUtils.hpp>
#include <polyfem/State.hpp>
#include "binding.hpp"
#include <pybind11/eigen.h>
#include <pybind11_json/pybind11_json.hpp>

namespace py = pybind11;
using namespace polyfem;
using namespace polyfem::mesh;
using namespace polyfem::solver;

void define_opt_utils(py::module_ &m)
{
  m.def(
       "apply_slim",
       [](const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
          const Eigen::MatrixXd &Vnew) {
         Eigen::MatrixXd Vsmooth;
         bool succeed = apply_slim(V, F, Vnew, Vsmooth, 1000);
         if (!succeed)
           throw std::runtime_error("SLIM failed to converge!");
         return Vsmooth;
       },
       py::arg("Vold"), py::arg("faces"), py::arg("Vnew"))

      .def("map_primitive_to_node_order",
           &AdjointTools::map_primitive_to_node_order, py::arg("state"),
           py::arg("primitives"))

      .def("map_node_to_primitive_order",
           &AdjointTools::map_node_to_primitive_order, py::arg("state"),
           py::arg("nodes"));
}
