#include <polyfem/State.hpp>
#include <polyfem/io/YamlToJson.hpp>
#include <polyfem/varforms/VarForm.hpp>

#include <algorithm>
#include <filesystem>
#include <fstream>
#include <stdexcept>
#include <string>
#include <vector>

#include "binding_wrapper.hpp"
#include "array_mesh_resources.hpp"

using namespace polyfem;

namespace
{

class PDEs
{
};

bool load_json(const std::string &json_file, json &out)
{
  std::ifstream file(json_file);
  if (!file.is_open())
    return false;

  file >> out;
  if (!out.contains("root_path"))
    out["root_path"] = json_file;
  return true;
}

bool load_yaml(const std::string &yaml_file, json &out)
{
  try
  {
    const io::FileSystemIO resources(std::filesystem::current_path());
    out = io::yaml_string_to_json(resources.read_string(yaml_file));
    if (!out.contains("root_path"))
      out["root_path"] = yaml_file;
  }
  catch (...)
  {
    return false;
  }
  return true;
}

void init_globals(State &state)
{
  static bool initialized = false;
  if (initialized)
    return;

  state.set_max_threads(1);
  state.init_logger(
      "",
      spdlog::level::level_enum::info,
      spdlog::level::level_enum::debug,
      false);
  initialized = true;
}

json settings_from_python(const py::object &settings)
{
  const std::string json_string = nb::cast<std::string>(py::str(settings));
  return json::parse(json_string);
}

class PythonSolver
{
public:
  PythonSolver() = default;

  void set_settings(const py::object &settings, bool strict_validation)
  {
    init_globals(state_);
    state_.init(settings_from_python(settings), strict_validation);
  }

  void set_settings_with_mesh(const py::object &settings,
                              const Eigen::MatrixXd &vertices,
                              const Eigen::MatrixXi &cells)
  {
    init_globals(state_);
    json args = settings_from_python(settings);
    json &geometry = args["geometry"];
    if (geometry.is_object()) geometry = json::array({geometry});
    json *body = nullptr;
    for (auto &entry : geometry)
    {
      if (!entry.value("enabled", true) || entry.value("is_obstacle", false)) continue;
      if (body != nullptr || entry.value("type", std::string("mesh")) != "mesh")
        throw std::invalid_argument("Array input requires one enabled deformable mesh.");
      body = &entry;
    }
    if (body == nullptr) throw std::invalid_argument("Array input requires one enabled deformable mesh.");
    (*body)["mesh"] = ArrayMeshResources::mesh_path();
    const ArrayMeshResources resources(vertices, cells);
    state_.init(args, resources, true);
    state_.load_mesh();
  }

  std::string settings() const
  {
    return state_.args.dump(2);
  }

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

  void set_mesh(
      const Eigen::MatrixXd &vertices,
      const Eigen::MatrixXi &cells,
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

void run_polyfem_command(
    const std::string &json_file,
    const std::string &yaml_file,
    const int log_level,
    const bool strict_validation,
    const int max_threads,
    const std::string &output_dir)
{
  json in_args = json::object();
  const bool ok = !json_file.empty()
                      ? load_json(json_file, in_args)
                      : load_yaml(yaml_file, in_args);

  if (!ok)
  {
    const std::string target = !json_file.empty() ? json_file : yaml_file;
    throw std::runtime_error("unable to open " + target + " file");
  }

  json patch = json::object();
  patch["/output/log/level"_json_pointer] = log_level;
  patch["/solver/max_threads"_json_pointer] = max_threads;
  if (!output_dir.empty())
    patch["/output/directory"_json_pointer] = std::filesystem::absolute(output_dir);
  in_args.merge_patch(patch);

  State state;
  init_globals(state);
  state.init(in_args, strict_validation);
  state.load_mesh();
  state.set_log_level(static_cast<spdlog::level::level_enum>(
      std::max(0, std::min(6, log_level))));

  Eigen::MatrixXd sol;
  {
    py::gil_scoped_release release;
    state.solve(sol);

    if (state.variational_formulation)
    {
      state.variational_formulation->compute_errors(sol);
      state.variational_formulation->save_json(sol);
      state.variational_formulation->export_data(sol);
    }
  }
}

} // namespace

void define_pde_types(py::module_ &m)
{
  const auto &pdes = py::class_<PDEs>(m, "PDEs");

  const std::vector<std::string> materials = {
      "LinearElasticity",
      "HookeLinearElasticity",
      "SaintVenant",
      "NeoHookean",
      "MooneyRivlin",
      "MooneyRivlin3Param",
      "MooneyRivlin3ParamSymbolic",
      "UnconstrainedOgden",
      "IncompressibleOgden",
      "Stokes",
      "NavierStokes",
      "OperatorSplitting",
      "IncompressibleLinearElasticity",
      "Laplacian",
      "Helmholtz",
      "Bilaplacian",
      "AMIPS",
      "FixedCorotational"};

  for (const auto &material : materials)
    pdes.attr(material.c_str()) = material;

  pdes.doc() = "List of supported partial differential equations";

  m.def(
      "is_tensor",
      [](const std::string &pde) {
        return !(pde == "Laplacian" || pde == "Helmholtz" || pde == "Bilaplacian");
      },
      "returns true if the pde is tensorial",
      py::arg("pde"));
}

void define_solver(py::module_ &m)
{
  py::class_<PythonSolver>(m, "Solver")
      .def(py::init<>())
      .def("set_settings_with_mesh", &PythonSolver::set_settings_with_mesh,
           "load settings and an array mesh through the geometry loader",
           py::arg("json"), py::arg("vertices"), py::arg("connectivity"))
      .def("settings", &PythonSolver::settings, "get solver settings as JSON text")
      .def(
          "set_settings",
          &PythonSolver::set_settings,
          "load solver settings from JSON text",
          py::arg("json"),
          py::arg("strict_validation") = false)
      .def(
          "set_max_threads",
          &PythonSolver::set_max_threads,
          "set maximum number of threads",
          py::arg("nthreads"))
      .def(
          "set_log_level",
          &PythonSolver::set_log_level,
          "set PolyFEM log level, from 0 (all logs) to 6 (no logs)",
          py::arg("log_level"))
      .def(
          "load_mesh_from_settings",
          &PythonSolver::load_mesh_from_settings,
          "load the mesh described by the current JSON settings")
      .def(
          "set_mesh",
          &PythonSolver::set_mesh,
          "load a mesh from vertices and connectivity",
          py::arg("vertices"),
          py::arg("connectivity"),
          py::arg("non_conforming") = false)
      .def(
          "solve",
          &PythonSolver::solve,
          "solve using the VarForm backend and return the raw solution bundle",
          py::arg("log_level") = 3);
}

void define_solve(py::module_ &m)
{
  m.def(
      "polyfem_command",
      &run_polyfem_command,
      "runs the PolyFEM command entry point",
      py::kw_only(),
      py::arg("json"),
      py::arg("yaml") = std::string(""),
      py::arg("log_level") = int(1),
      py::arg("strict_validation") = bool(true),
      py::arg("max_threads") = int(1),
      py::arg("output_dir") = "");
}
