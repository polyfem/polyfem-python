#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <stdexcept>

#include <pybind11_json/pybind11_json.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include "differentiable/binding.hpp"
#include "mesh/binding.hpp"
#include "state/binding.hpp"
#include "solver/binding.hpp"

namespace py = pybind11;

PYBIND11_MODULE(polyfempy, m)
{
  define_pde_types(m);

  define_solver(m);
  define_solve(m);

  define_mesh(m);

  define_nonlinear_problem(m);

  define_differentiable_cache(m);
  define_adjoint(m);
  define_objective(m);
  define_opt_utils(m);
}
