#include <polyfem/solver/NLProblem.hpp>
#include "binding.hpp"
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace polyfem;
using namespace polyfem::solver;

void define_nonlinear_problem(py::module_ &m)
{
  py::class_<FullNLProblem>(m, "FullNLProblem",
                            "Full nonlinear problem in the simulation")
      .def("init", &FullNLProblem::init, "Initialization", py::arg("x0"))
      .def("value", &FullNLProblem::value)
      .def(
          "gradient",
          [](FullNLProblem &prob, const Eigen::VectorXd &x) {
            Eigen::VectorXd grad;
            prob.gradient(x, grad);
            return grad;
          },
          py::arg("x"))
      .def(
          "hessian",
          [](FullNLProblem &prob, const Eigen::VectorXd &x) {
            StiffnessMatrix hess;
            prob.hessian(x, hess);
            return hess;
          },
          py::arg("x"))
      .def("is_step_valid", &FullNLProblem::is_step_valid, py::arg("x0"),
           py::arg("x1"))
      .def("is_step_collision_free", &FullNLProblem::is_step_collision_free,
           py::arg("x0"), py::arg("x1"))
      .def("max_step_size", &FullNLProblem::max_step_size, py::arg("x0"),
           py::arg("x1"))
      .def("line_search_begin", &FullNLProblem::line_search_begin,
           py::arg("x0"), py::arg("x1"))
      .def("line_search_end", &FullNLProblem::line_search_end)
      .def("solution_changed", &FullNLProblem::solution_changed, py::arg("x"))
      .def("stop", &FullNLProblem::stop, py::arg("x"));

  py::class_<NLProblem, FullNLProblem>(m, "NLProblem", "Nonlinear problem in the simulation")
  .def("full_to_reduced", &NLProblem::full_to_reduced, py::arg("full"))
  .def("reduced_to_full", &NLProblem::reduced_to_full, py::arg("reduced"));
}