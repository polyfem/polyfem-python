// #include <polyfem/solver/AdjointTools.hpp>
#include <polyfem/solver/forms/adjoint_forms/AdjointForm.hpp>
#include <polyfem/solver/forms/adjoint_forms/VariableToSimulation.hpp>
#include <polyfem/State.hpp>
#include <polyfem/solver/Optimizations.hpp>
#include <polyfem/utils/MatrixUtils.hpp>
#include <polyfem/State.hpp>
#include "binding.hpp"
#include <pybind11/eigen.h>
#include <pybind11_json/pybind11_json.hpp>

namespace py = pybind11;
using namespace polyfem;
using namespace polyfem::solver;

void define_objective(py::module_ &m)
{
  py::class_<AdjointForm, std::shared_ptr<AdjointForm>>(m, "Objective")
      .def("name", &AdjointForm::name)

      .def("value", &AdjointForm::value, py::arg("x"))

      .def("solution_changed", &AdjointForm::solution_changed, py::arg("x"))

      .def("is_step_collision_free", &AdjointForm::is_step_collision_free,
           py::arg("x0"), py::arg("x1"))

      .def("max_step_size", &AdjointForm::max_step_size, py::arg("x0"),
           py::arg("x1"))

      .def(
          "derivative",
          [](AdjointForm &obj, State &solver, const Eigen::VectorXd &x,
             const std::string &wrt) -> Eigen::VectorXd {
            if (wrt == "solution")
              return obj.compute_adjoint_rhs(x, solver);
            else if (wrt == obj.get_variable_to_simulations()[0]->name())
            {
              Eigen::VectorXd grad;
              obj.compute_partial_gradient(x, grad);
              return grad;
            }
            else
              throw std::runtime_error(
                  "Input type does not match objective derivative type!");
          },
          py::arg("solver"), py::arg("x"), py::arg("wrt"));

  m.def("create_objective", &AdjointOptUtils::create_simple_form,
        py::arg("obj_type"), py::arg("param_type"), py::arg("solver"),
        py::arg("parameters"));
}
