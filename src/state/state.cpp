#include <polyfem/State.hpp>
#include <polyfem/varforms/VarForm.hpp>

#include <algorithm>
#include <stdexcept>
#include <string>

#include "binding_wrapper.hpp"
#include "array_mesh_resources.hpp"

using namespace polyfem;

namespace
{

  void init_globals(State &state)
  {
    static bool initialized = false;
    if (initialized)
      return;

    state.set_max_threads(1);
    state.init_logger("", spdlog::level::level_enum::info,
                      spdlog::level::level_enum::debug, false);
    initialized = true;
  }

  class PythonSolver
  {
public:
    PythonSolver() = default;

    void set_settings(const std::string &settings, bool strict_validation)
    {
      init_globals(state_);
      state_.init(json::parse(settings), strict_validation);
    }

    void set_settings_with_mesh(const std::string &settings,
                                const Eigen::MatrixXd &vertices,
                                const Eigen::MatrixXi &cells)
    {
      init_globals(state_);
      json args = json::parse(settings);
      json &geometry = args["geometry"];
      if (geometry.is_object())
        geometry = json::array({geometry});
      json *body = nullptr;
      for (auto &entry : geometry)
      {
        if (!entry.value("enabled", true) || entry.value("is_obstacle", false))
          continue;
        if (body != nullptr
            || entry.value("type", std::string("mesh")) != "mesh")
          throw std::invalid_argument(
              "Array input requires one enabled deformable mesh.");
        body = &entry;
      }
      if (body == nullptr)
        throw std::invalid_argument(
            "Array input requires one enabled deformable mesh.");
      (*body)["mesh"] = ArrayMeshResources::mesh_path();
      const ArrayMeshResources resources(vertices, cells);
      state_.init(args, resources, true);
      state_.load_mesh();
    }

    std::string settings() const { return state_.args.dump(2); }

    void set_max_threads(const int nthreads)
    {
      state_.set_max_threads(nthreads);
    }

    void set_log_level(int log_level)
    {
      init_globals(state_);
      log_level = std::max(0, std::min(6, log_level));
      state_.set_log_level(static_cast<spdlog::level::level_enum>(log_level));
    }

    void load_mesh_from_settings()
    {
      init_globals(state_);
      state_.load_mesh();
    }

    void set_mesh(const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &cells,
                  const bool non_conforming)
    {
      init_globals(state_);
      state_.set_mesh(vertices, cells, non_conforming);
    }

    py::dict solve(int log_level)
    {
      set_log_level(log_level);

      Eigen::MatrixXd sol;
      {
        // Worker threads must be able to call Python expressions.
        py::gil_scoped_release release;
        state_.solve(sol);

        if (state_.variational_formulation)
        {
          state_.variational_formulation->compute_errors(sol);
          state_.variational_formulation->save_json(sol);
          state_.variational_formulation->export_data(sol);
        }
      }

      py::dict meta;
      meta["from_bundle"] = true;
      meta["solution"] = "raw_backend_solution";

      py::dict bundle;
      bundle["_result_bundle"] = true;
      bundle["sol"] = sol;
      bundle["meta"] = meta;
      return bundle;
    }

private:
    State state_;
  };

} // namespace

void define_solver(py::module_ &m)
{
  py::class_<PythonSolver>(m, "Solver")
      .def(py::init<>())
      .def("set_settings_with_mesh", &PythonSolver::set_settings_with_mesh,
           "load settings and an array mesh through the geometry loader",
           py::arg("json"), py::arg("vertices"), py::arg("connectivity"))
      .def("settings", &PythonSolver::settings,
           "get solver settings as JSON text")
      .def("set_settings", &PythonSolver::set_settings,
           "load solver settings from JSON text", py::arg("json"),
           py::arg("strict_validation") = false)
      .def("set_max_threads", &PythonSolver::set_max_threads,
           "set maximum number of threads", py::arg("nthreads"))
      .def("set_log_level", &PythonSolver::set_log_level,
           "set PolyFEM log level, from 0 (all logs) to 6 (no logs)",
           py::arg("log_level"))
      .def("load_mesh_from_settings", &PythonSolver::load_mesh_from_settings,
           "load the mesh described by the current JSON settings")
      .def("set_mesh", &PythonSolver::set_mesh,
           "load a mesh from vertices and connectivity", py::arg("vertices"),
           py::arg("connectivity"), py::arg("non_conforming") = false)
      .def("solve", &PythonSolver::solve,
           "solve using the VarForm backend and return the raw solution bundle",
           py::arg("log_level") = 3);
}
