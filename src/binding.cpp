#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <stdexcept>

#include "binding_wrapper.hpp"

#include "mesh/binding.hpp"
#include "state/binding.hpp"

PY_MODULE(polyfempy, m)
{
  define_pde_types(m);

  define_solver(m);
  define_solve(m);

  define_mesh(m);

  m.def("version", []() { return "polyfempy nanobind backend"; }, "Get version information");
}
