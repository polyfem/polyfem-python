#include <polyfem/State.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/varforms/VarForm.hpp>

#include <nanobind/nanobind.h>
#include <nanobind/eigen/dense.h>
#include <nanobind/stl/string.h>

namespace py = nanobind;
using namespace polyfem;
using namespace py::literals;

NB_MODULE(polyfempy, m)
{
    py::class_<State>(m, "Solver")
        .def(py::init<>())
        .def("set_settings", [](State &state, const std::string &settings, bool strict) {
            state.init(json::parse(settings), strict);
        }, "json"_a, "strict_validation"_a = false)
        .def("set_settings_with_mesh", [](State &state, const std::string &settings,
                                         const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &cells) {
            state.init(json::parse(settings), true);
            const mesh::MeshData mesh(vertices, cells);
            state.load_mesh(false, &mesh);
        }, "json"_a, "vertices"_a, "connectivity"_a)
        .def("load_mesh_from_settings", [](State &state) { state.load_mesh(); })
        .def("set_mesh", &State::set_mesh,
             "vertices"_a, "connectivity"_a, "non_conforming"_a = false)
        .def("solve", [](State &state) {
            Eigen::MatrixXd sol;
            {
                // Python expressions may run on worker threads.
                py::gil_scoped_release release;
                state.solve(sol);
                state.variational_formulation->compute_errors(sol);
                const auto timings = state.variational_formulation->output_timings();
                logger().info("total time: {}s", timings.total_time());
                state.variational_formulation->save_json(sol);
                state.variational_formulation->export_data(sol);
            }
            return sol;
        });
}
