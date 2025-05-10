#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_nonlinear_problem(py::module_ &m);
