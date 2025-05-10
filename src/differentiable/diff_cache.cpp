#include <polyfem/solver/DiffCache.hpp>
#include <polyfem/State.hpp>
#include "binding.hpp"
#include <pybind11/eigen.h>

namespace py = pybind11;
using namespace polyfem::solver;

void define_differentiable_cache(py::module_ &m)
{
  py::enum_<CacheLevel>(m, "CacheLevel", "Caching level of the simulator.")
      .value("None", CacheLevel::None, "No cache at all")
      .value("Solution", CacheLevel::Solution, "Cache solutions")
      .value("Derivatives", CacheLevel::Derivatives,
             "Cache solutions and quantities for gradient computation")
      .export_values();

  py::class_<DiffCache>(m, "DiffCache", "Cache of the simulator")

      .def("size", &DiffCache::size,
           "Get current cache size (number of time steps)")

      .def("solution", &DiffCache::u, "Get solution",
           py::arg("time_step") = int(0))

      .def("displacement", &DiffCache::u, "Get displacement",
           py::arg("time_step") = int(0))

      .def("velocity", &DiffCache::v, "Get velocity",
           py::arg("time_step") = int(0))

      .def("acceleration", &DiffCache::acc, "Get acceleration",
           py::arg("time_step") = int(0))

      .def("hessian", &DiffCache::gradu_h, "Get energy hessian at solution",
           py::arg("time_step") = int(0));
}
