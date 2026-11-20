// nanobind binding wrapper
// pybind11 support has been removed - only nanobind is used

#pragma once

// Hard fail guard: detect if pybind11 was included before this header
#ifdef PYBIND11_VERSION
    #error "pybind11 headers included. pybind11 support has been removed - only nanobind is used."
#endif

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/function.h>
#include <nanobind/stl/vector.h>
#include <nanobind/stl/string.h>
#include <nanobind/stl/optional.h>
#include <nanobind/stl/shared_ptr.h>
#include <nanobind/ndarray.h>

namespace py = nanobind;
namespace nb = nanobind;
#define PY_MODULE(name, m) NB_MODULE(name, m)

// Compatibility macros for return_value_policy
// Note: Use 'inline constexpr' (not 'static inline constexpr') for namespace-scope variables in C++17
namespace py_return_value_policy {
    inline constexpr nanobind::rv_policy reference = nanobind::rv_policy::reference_internal;
    inline constexpr nanobind::rv_policy copy = nanobind::rv_policy::copy;
    inline constexpr nanobind::rv_policy move = nanobind::rv_policy::move;
    inline constexpr nanobind::rv_policy reference_internal = nanobind::rv_policy::reference_internal;
}

