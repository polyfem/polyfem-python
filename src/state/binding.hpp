#pragma once

#include "binding_wrapper.hpp"

void define_pde_types(py::module_ &m);
void define_solver(py::module_ &m);
void define_solve(py::module_ &m);
