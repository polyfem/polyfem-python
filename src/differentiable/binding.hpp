#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_differentiable_cache(py::module_ &m);
void define_adjoint(py::module_ &m);
void define_objective(py::module_ &m);
void define_opt_utils(py::module_ &m);
