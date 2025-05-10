#pragma once

#include <pybind11/pybind11.h>

namespace py = pybind11;

void define_pde_types(py::module_ &m);
void define_solver(py::module_ &m);
void define_solve(py::module_ &m);
