#include <polyfem/mesh/MeshUtils.hpp>
#include <polyfem/assembler/AssemblerUtils.hpp>
#include <polyfem/assembler/GenericProblem.hpp>
#include <polyfem/io/Evaluator.hpp>
#include <polyfem/io/YamlToJson.hpp>
#include <polyfem/utils/Logger.hpp>
#include <polyfem/utils/StringUtils.hpp>
#include <polyfem/utils/JSONUtils.hpp>
#include <polyfem/utils/GeogramUtils.hpp>
#include <polyfem/solver/NLProblem.hpp>
#include <polyfem/time_integrator/ImplicitTimeIntegrator.hpp>
#include <polyfem/State.hpp>

#include <stdexcept>

#include <pybind11_json/pybind11_json.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

#include "differentiable/binding.hpp"

namespace py = pybind11;
using namespace polyfem;

typedef std::function<Eigen::MatrixXd(double x, double y, double z)> BCFuncV;
typedef std::function<double(double x, double y, double z)> BCFuncS;

class Assemblers
{
};

class PDEs
{
};

// TODO add save_time_sequence

namespace
{

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
      out = io::yaml_file_to_json(yaml_file);
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

    if (!initialized)
    {
      state.set_max_threads(1);
      state.init_logger("", spdlog::level::level_enum::info,
                        spdlog::level::level_enum::debug, false);

      initialized = true;
    }
  }

} // namespace

void define_pde_types(py::module_ &m)
{
  const auto &pdes = py::class_<PDEs>(m, "PDEs");

  const std::vector<std::string> materials = {"LinearElasticity",
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

  for (const auto &a : materials)
    pdes.attr(a.c_str()) = a;

  pdes.doc() = "List of supported partial differential equations";

  m.def(
      "is_tensor",
      [](const std::string &pde) {
        if (pde == "Laplacian" || pde == "Helmholtz" || pde == "Bilaplacian")
          return false;
        return true;
      },
      "returns true if the pde is tensorial", py::arg("pde"));
}

void define_solver(py::module_ &m)
{
  const auto setting_lambda = [](State &self, const py::object &settings,
                                 bool strict_validation) {
    using namespace polyfem;

    init_globals(self);
    // py::scoped_ostream_redirect output;
    const std::string json_string = py::str(settings);
    self.init(json::parse(json_string), strict_validation);
  };

  py::class_<State, std::shared_ptr<State>>(m, "Solver")
      .def(py::init<>())

      .def("is_tensor", [](const State &s) { return s.assembler->is_tensor(); })

      .def(
          "settings", [](const State &s) { return s.args; },
          "get PDE and problem parameters from the solver")

      .def("set_settings", setting_lambda,
           "load PDE and problem parameters from the settings", py::arg("json"),
           py::arg("strict_validation") = false)

      .def("set_max_threads", &State::set_max_threads,
           "set maximum number of threads", py::arg("nthreads"))

      .def("ndof", &State::ndof, "Dimension of the solution")

      .def(
          "n_bases", [](const State &s) { return s.n_bases; },
          "Number of basis")

      .def(
          "set_log_level",
          [](State &s, int log_level) {
            init_globals(s);
            //    py::scoped_ostream_redirect output;
            log_level = std::max(0, std::min(6, log_level));
            s.set_log_level(static_cast<spdlog::level::level_enum>(log_level));
          },
          "sets polyfem log level, valid value between 0 (all logs) and 6 (no logs)",
          py::arg("log_level"))

      .def(
          "mesh", [](State &s) -> mesh::Mesh & { return *s.mesh.get(); },
          "Get mesh in simulator", py::return_value_policy::reference)

      .def(
          "load_mesh_from_settings",
          [](State &s) {
            init_globals(s);
            s.load_mesh();
          },
          "Loads a mesh from the 'mesh' field of the json")

      // .def(
      //     "reload_boundary_conditions",
      //     [](State &s) {
      //       auto bc = s.args["boundary_conditions"];
      //       bc["root_path"] = s.root_path();
      //       s.problem->clear();
      //       s.problem->set_parameters(bc);
      //     },
      //     "Reload boundary conditions from the json.")

      .def(
          "update_dirichlet_nodes",
          [](State &s, const Eigen::VectorXi &node_ids,
             const Eigen::MatrixXd &nodal_dirichlet) {
            auto tensor_problem = std::dynamic_pointer_cast<
                polyfem::assembler::GenericTensorProblem>(s.problem);
            // if (!s.iso_parametric())
            //   throw std::runtime_error(
            //       "Can only update dirichlet nodes for isoparametric.");
            tensor_problem->update_dirichlet_nodes(s.in_node_to_node, node_ids,
                                                   nodal_dirichlet);
          },
          "Reload boundary conditions from the json.")

      .def(
          "load_mesh_from_path",
          [](State &s, const std::string &path, const bool normalize_mesh,
             const double vismesh_rel_area, const int n_refs,
             const double boundary_id_threshold) {
            init_globals(s);
            s.args["geometry"] = R"([{ }])"_json;
            s.args["geometry"][0]["mesh"] = path;
            s.args["geometry"][0]["advanced"]["normalize_mesh"] =
                normalize_mesh;
            s.args["geometry"][0]["surface_selection"] =
                R"({ "threshold": 0.0 })"_json;
            s.args["geometry"][0]["surface_selection"]["threshold"] =
                boundary_id_threshold;
            s.args["geometry"][0]["n_refs"] = n_refs;
            s.args["output"]["paraview"]["vismesh_rel_area"] = vismesh_rel_area;
            s.load_mesh();
          },
          "Loads a mesh from the path and 'bc_tag' from the json if any bc tags",
          py::arg("path"), py::arg("normalize_mesh") = bool(false),
          py::arg("vismesh_rel_area") = double(0.00001),
          py::arg("n_refs") = int(0),
          py::arg("boundary_id_threshold") = double(-1))

      .def(
          "load_mesh_from_path_and_tags",
          [](State &s, const std::string &path, const std::string &bc_tag,
             const bool normalize_mesh, const double vismesh_rel_area,
             const int n_refs, const double boundary_id_threshold) {
            init_globals(s);
            s.args["geometry"] = R"([{ }])"_json;
            s.args["geometry"][0]["mesh"] = path;
            s.args["bc_tag"] = bc_tag;
            s.args["geometry"][0]["advanced"]["normalize_mesh"] =
                normalize_mesh;
            s.args["geometry"][0]["surface_selection"] =
                R"({ "threshold": 0.0 })"_json;
            s.args["geometry"][0]["surface_selection"]["threshold"] =
                boundary_id_threshold;
            s.args["geometry"][0]["n_refs"] = n_refs;
            s.args["output"]["paraview"]["vismesh_rel_area"] = vismesh_rel_area;
            s.load_mesh();
          },
          "Loads a mesh and bc_tags from path", py::arg("path"),
          py::arg("bc_tag_path"), py::arg("normalize_mesh") = bool(false),
          py::arg("vismesh_rel_area") = double(0.00001),
          py::arg("n_refs") = int(0),
          py::arg("boundary_id_threshold") = double(-1))

      .def(
          "set_mesh",
          [](State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
             const int n_refs, const double boundary_id_threshold) {
            init_globals(s);
            s.mesh = mesh::Mesh::create(V, F);
            s.args["geometry"] = R"([{ }])"_json;
            s.args["geometry"][0]["n_refs"] = n_refs;
            s.args["geometry"][0]["surface_selection"] =
                R"({ "threshold": 0.0 })"_json;
            s.args["geometry"][0]["surface_selection"]["threshold"] =
                boundary_id_threshold;

            s.load_mesh();
          },
          "Loads a mesh from vertices and connectivity", py::arg("vertices"),
          py::arg("connectivity"), py::arg("n_refs") = int(0),
          py::arg("boundary_id_threshold") = double(-1))

      // .def(
      //     "set_mesh",
      //     [](State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
      //        const std::string &surface_selections_file,
      //        const int volume_selection) {
      //       init_globals(s);
      //       s.mesh = mesh::Mesh::create(V, F);
      //       s.args["geometry"] = R"([{ }])"_json;
      //       s.args["geometry"][0]["surface_selection"] =
      //           surface_selections_file;
      //       s.args["geometry"][0]["volume_selection"] = volume_selection;

      //       s.load_mesh();
      //     },
      //     "Loads a mesh from vertices and connectivity, specifying surfaces",
      //     py::arg("vertices"), py::arg("connectivity"),
      //     py::arg("surface_selections_file") = "",
      //     py::arg("volume_selection") = int(1))

      .def(
          "set_high_order_mesh",
          [](State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F,
             const Eigen::MatrixXd &nodes_pos,
             const std::vector<std::vector<int>> &nodes_indices,
             const bool normalize_mesh, const double vismesh_rel_area,
             const int n_refs, const double boundary_id_threshold) {
            init_globals(s);
            //    py::scoped_ostream_redirect output;

            s.mesh = mesh::Mesh::create(V, F);
            s.mesh->attach_higher_order_nodes(nodes_pos, nodes_indices);

            s.args["geometry"][0]["advanced"]["normalize_mesh"] =
                normalize_mesh;
            s.args["geometry"][0]["n_refs"] = n_refs;
            s.args["geometry"][0]["surface_selection"] =
                R"({ "threshold": 0.0 })"_json;
            s.args["geometry"][0]["surface_selection"]["threshold"] =
                boundary_id_threshold;
            s.args["output"]["paraview"]["vismesh_rel_area"] = vismesh_rel_area;

            s.load_mesh();
          },
          "Loads an high order mesh from vertices, connectivity, nodes, and node indices mapping element to nodes",
          py::arg("vertices"), py::arg("connectivity"), py::arg("nodes_pos"),
          py::arg("nodes_indices"), py::arg("normalize_mesh") = bool(false),
          py::arg("vismesh_rel_area") = double(0.00001),
          py::arg("n_refs") = int(0),
          py::arg("boundary_id_threshold") = double(-1))

      .def(
          "get_vertices",
          [](const State &state) {
            Eigen::MatrixXd vertices;
            state.get_vertices(vertices);
            return vertices;
          },
          "get the vertices")

      .def(
          "get_elements",
          [](const State &state) {
            Eigen::MatrixXi elements;
            state.get_elements(elements);
            return elements;
          },
          "get the elements")

      .def(
          "nl_problem", [](State &s) { return *(s.solve_data.nl_problem); },
          py::return_value_policy::reference)

      .def(
          "solve",
          [](State &s, int log_level) {
            init_globals(s);
            //    py::scoped_ostream_redirect output;
            s.stats.compute_mesh_stats(*s.mesh);

            s.build_basis();

            s.assemble_rhs();
            s.assemble_mass_mat();

            s.set_log_level(static_cast<spdlog::level::level_enum>(log_level));

            Eigen::MatrixXd sol, pressure;
            s.solve_problem(sol, pressure);

            s.compute_errors(sol);

            s.save_json(sol);
            s.export_data(sol, pressure);

            return py::make_tuple(sol, pressure);
          },
          "solve the pde", py::arg("log_level") = int(3))
      .def(
          "build_basis",
          [](State &s) {
            if (!s.mesh)
              throw std::runtime_error("Load mesh first!");

            s.build_basis();
          },
          "build finite element basis")
      .def(
          "assemble",
          [](State &s) {
            if (s.bases.size() == 0)
              throw std::runtime_error("Call build_basis() first!");

            s.assemble_rhs();
            s.assemble_mass_mat();
          },
          "assemble RHS and mass matrix if needed")
      .def(
          "init_timestepping",
          [](State &s, const double t0, const double dt) {
            if (!s.solve_data.rhs_assembler || s.mass.size() == 0)
              throw std::runtime_error("Call assemble() first!");

            s.solution_frames.clear();
            Eigen::MatrixXd sol, pressure;
            s.init_solve(sol, pressure);
            s.init_nonlinear_tensor_solve(sol, t0 + dt);
            s.cache_transient_adjoint_quantities(
                0, sol,
                Eigen::MatrixXd::Zero(s.mesh->dimension(),
                                      s.mesh->dimension()));
            return sol;
          },
          "initialize timestepping", py::arg("t0"), py::arg("dt"))
      .def(
          "step_in_time",
          [](State &s, Eigen::MatrixXd &sol, const double t0, const double dt,
             const int t) {
            if (s.assembler->name() == "NavierStokes"
                || s.assembler->name() == "OperatorSplitting"
                || s.is_problem_linear() || s.is_homogenization())
              throw std::runtime_error("Formulation " + s.assembler->name()
                                       + " is not supported!");

            s.solve_tensor_nonlinear(sol, t);
            s.cache_transient_adjoint_quantities(
                t, sol,
                Eigen::MatrixXd::Zero(s.mesh->dimension(),
                                      s.mesh->dimension()));

            s.solve_data.time_integrator->update_quantities(sol);
            s.solve_data.nl_problem->update_quantities(t0 + (t + 1) * dt, sol);
            s.solve_data.update_dt();
            s.solve_data.update_barrier_stiffness(sol);
            return sol;
          },
          "step in time", py::arg("solution"), py::arg("t0"), py::arg("dt"),
          py::arg("t"))

      .def(
          "solve_adjoint",
          [](State &s, const Eigen::MatrixXd &adjoint_rhs) {
            if (adjoint_rhs.cols() != s.diff_cached.size()
                || adjoint_rhs.rows() != s.diff_cached.u(0).size())
              throw std::runtime_error("Invalid adjoint_rhs shape!");
            if (!s.problem->is_time_dependent() && !s.lin_solver_cached
                && s.is_homogenization()) // nonlinear static solve only
            {
              Eigen::MatrixXd reduced;
              for (int i = 0; i < adjoint_rhs.cols(); i++)
              {
                Eigen::VectorXd reduced_vec =
                    s.solve_data.nl_problem->full_to_reduced_grad(
                        adjoint_rhs.col(i));
                if (i == 0)
                  reduced.setZero(reduced_vec.rows(), adjoint_rhs.cols());
                reduced.col(i) = reduced_vec;
              }
              return s.solve_adjoint_cached(reduced);
            }
            else
              return s.solve_adjoint_cached(adjoint_rhs);
          },
          "Solve the adjoint equation given the gradient of objective wrt. PDE solution")

      .def(
          "set_cache_level",
          [](State &s, solver::CacheLevel level) {
            s.optimization_enabled = level;
            if (level == solver::CacheLevel::Derivatives)
            {
              if (s.is_contact_enabled())
              {
                if (!s.args["contact"]["use_convergent_formulation"])
                {
                  s.args["contact"]["use_convergent_formulation"] = true;
                  logger().info(
                      "Use convergent formulation for differentiable contact...");
                }
                if (s.args["/solver/contact/barrier_stiffness"_json_pointer]
                        .is_string())
                {
                  logger().error(
                      "Only constant barrier stiffness is supported in differentiable contact!");
                }
              }

              if (s.args.contains("boundary_conditions")
                  && s.args["boundary_conditions"].contains("rhs"))
              {
                json rhs = s.args["boundary_conditions"]["rhs"];
                if ((rhs.is_array() && rhs.size() > 0 && rhs[0].is_string())
                    || rhs.is_string())
                  logger().error(
                      "Only constant rhs over space is supported in differentiable code!");
              }
            }
          },
          "Set solution caching level", py::arg("cache_level"))

      .def(
          "get_solution_cache", [](State &s) { return s.diff_cached; },
          "get the cached solution after simulation, this function requires setting CacheLevel before the simulation")

      .def("get_solutions",
           [](State &s) {
             Eigen::MatrixXd sol(s.diff_cached.u(0).size(),
                                 s.diff_cached.size());
             for (int i = 0; i < sol.cols(); i++)
               sol.col(i) = s.diff_cached.u(i);
             return sol;
           })

      .def(
          "compute_errors",
          [](State &s, Eigen::MatrixXd &sol) { s.compute_errors(sol); },
          "compute the error", py::arg("solution"))

      .def(
          "export_data",
          [](State &s, const Eigen::MatrixXd &sol,
             const Eigen::MatrixXd &pressure) { s.export_data(sol, pressure); },
          "exports all data specified in the settings")
      .def(
          "export_vtu",
          [](State &s, std::string &path, const Eigen::MatrixXd &sol,
             const Eigen::MatrixXd &pressure, const double time,
             const double dt) {
            s.out_geom.save_vtu(
                s.resolve_output_path(path), s, sol, pressure, time, dt,
                io::OutGeometryData::ExportOptions(s.args, s.mesh->is_linear(),
                                                   s.problem->is_scalar(),
                                                   s.solve_export_to_file),
                s.is_contact_enabled(), s.solution_frames);
          },
          "exports the solution as vtu", py::arg("path"), py::arg("solution"),
          py::arg("pressure") = Eigen::MatrixXd(), py::arg("time") = double(0.),
          py::arg("dt") = double(0.))
      .def(
          "set_friction_coefficient",
          [](State &self, const double mu) {
            self.args["contact"]["friction_coefficient"] = mu;
          },
          "set friction coefficient", py::arg("mu"))
      .def(
          "set_initial_velocity",
          [](State &self, const int body_id, const Eigen::VectorXd &velocity) {
            if (self.bases.size() == 0)
              log_and_throw_adjoint_error("Build basis first!");

            if (velocity.size() != self.mesh->dimension())
              log_and_throw_adjoint_error("Invalid velocity size {}!",
                                          velocity.size());

            // Initialize initial velocity
            if (self.initial_vel_update.size() != self.ndof())
              log_and_throw_adjoint_error("Call init_timestepping first!");

            assert(self.initial_vel_update.size() == self.ndof());
            // Set initial velocity
            for (size_t e = 0; e < self.bases.size(); e++)
            {
              if (self.mesh->get_body_id(e) == body_id)
              {
                const auto &bs = self.bases[e];
                for (const auto &b : bs.bases)
                  for (const auto &g : b.global())
                    for (int d = 0; d < velocity.size(); d++)
                      self.initial_vel_update(g.index * velocity.size() + d) =
                          velocity(d);
              }
            }
          },
          "set initial velocity for one body", py::arg("body_id"),
          py::arg("velocity"))
      .def(
          "set_initial_displacement",
          [](State &self, const int body_id, const Eigen::VectorXd &disp) {
            if (self.bases.size() == 0)
              log_and_throw_adjoint_error("Build basis first!");

            if (disp.size() != self.mesh->dimension())
              log_and_throw_adjoint_error("Invalid disp size {}!", disp.size());

            // Initialize initial displacement
            if (self.initial_sol_update.size() != self.ndof())
              log_and_throw_adjoint_error("Call init_timestepping first!");

            assert(self.initial_sol_update.size() == self.ndof());
            // Set initial displacement
            for (size_t e = 0; e < self.bases.size(); e++)
            {
              if (self.mesh->get_body_id(e) == body_id)
              {
                const auto &bs = self.bases[e];
                for (const auto &b : bs.bases)
                  for (const auto &g : b.global())
                    for (int d = 0; d < disp.size(); d++)
                      self.initial_sol_update(g.index * disp.size() + d) =
                          disp(d);
              }
            }
          },
          "set initial displacement for one body", py::arg("body_id"),
          py::arg("displacement"))
      .def(
          "set_per_element_material",
          [](State &self, const Eigen::VectorXd &lambda,
             const Eigen::VectorXd &mu) {
            if (self.bases.size() == 0)
              log_and_throw_adjoint_error("Build basis first!");

            assert(lambda.size() == self.bases.size());
            assert(mu.size() == self.bases.size());
            self.assembler->update_lame_params(lambda, mu);
          },
          "set per-element Lame parameters", py::arg("lambda"), py::arg("mu"));
}

void define_solve(py::module_ &m)
{

  m.def(
      "polyfem_command",
      [](const std::string &json_file, const std::string &yaml_file,
         const int log_level, const bool strict_validation,
         const int max_threads, const std::string &output_dir) {
        json in_args = json({});

        const bool ok = !json_file.empty() ? load_json(json_file, in_args)
                                           : load_yaml(yaml_file, in_args);

        if (!ok)
          throw std::runtime_error(
              fmt::format("unable to open {} file", json_file));

        json tmp = json::object();
        tmp["/output/log/level"_json_pointer] = int(log_level);
        tmp["/solver/max_threads"_json_pointer] = max_threads;
        if (!output_dir.empty())
          tmp["/output/directory"_json_pointer] =
              std::filesystem::absolute(output_dir);
        assert(tmp.is_object());
        in_args.merge_patch(tmp);

        std::vector<std::string> names;
        std::vector<Eigen::MatrixXi> cells;
        std::vector<Eigen::MatrixXd> vertices;

        State state;
        state.init(in_args, strict_validation);
        state.load_mesh(/*non_conforming=*/false, names, cells, vertices);

        // Mesh was not loaded successfully; load_mesh() logged the error.
        if (state.mesh == nullptr)
          throw std::runtime_error("Failed to load the mesh!");

        state.stats.compute_mesh_stats(*state.mesh);

        state.build_basis();

        state.assemble_rhs();
        state.assemble_mass_mat();

        Eigen::MatrixXd sol;
        Eigen::MatrixXd pressure;

        state.solve_problem(sol, pressure);

        state.compute_errors(sol);

        state.save_json(sol);
        state.export_data(sol, pressure);
      },
      "runs the polyfem command, internal usage", py::kw_only(),
      py::arg("json"), py::arg("yaml") = std::string(""),
      py::arg("log_level") = int(1), py::arg("strict_validation") = bool(true),
      py::arg("max_threads") = int(1), py::arg("output_dir") = "");

  //   m.def(
  //       "solve_febio",
  //       [](const std::string &febio_file, const std::string &output_path,
  //          const int log_level, const py::kwargs &kwargs) {
  //         if (febio_file.empty())
  //           throw pybind11::value_error("Specify a febio file!");

  //         // json in_args = opts.is_none() ? json({}) : json(opts);
  //         json in_args = json(static_cast<py::dict>(kwargs));

  //         if (!output_path.empty())
  //         {
  //           in_args["export"]["paraview"] = output_path;
  //           in_args["export"]["wire_mesh"] =
  //               utils::StringUtils::replace_ext(output_path, "obj");
  //           in_args["export"]["material_params"] = true;
  //           in_args["export"]["body_ids"] = true;
  //           in_args["export"]["contact_forces"] = true;
  //           in_args["export"]["surface"] = true;
  //         }

  //         const int discr_order =
  //             in_args.contains("discr_order") ? int(in_args["discr_order"]) :
  //             1;

  //         if (discr_order == 1 && !in_args.contains("vismesh_rel_area"))
  //           in_args["output"]["paraview"]["vismesh_rel_area"] = 1e10;

  //         State state;
  //         state.init_logger("", log_level, false);
  //         state.init(in_args);
  //         state.load_febio(febio_file, in_args);
  //         state.stats.compute_mesh_stats(*state.mesh);

  //         state.build_basis();

  //         state.assemble_rhs();
  //         state.assemble_mass_mat();

  //         Eigen::MatrixXd sol, pressure;
  //         state.solve_problem(sol, pressure);

  //         state.save_json();
  //         state.export_data(sol, pressure);
  //       },
  //       "runs FEBio", py::arg("febio_file"),
  //       py::arg("output_path") = std::string(""), py::arg("log_level") = 2);

  //   m.def(
  //       "solve",
  //       [](const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &cells,
  //          const py::object &sidesets_func, const int log_level,
  //          const py::kwargs &kwargs) {
  //         std::string log_file = "";

  //         std::unique_ptr<State> res =
  //             std::make_unique<State>();
  //         State &state = *res;
  //         state.init_logger(log_file, log_level, false);

  //         json in_args = json(static_cast<py::dict>(kwargs));

  //         state.init(in_args);

  //         state.load_mesh(vertices, cells);

  //         [&]() {
  //           if (!sidesets_func.is_none())
  //           {
  //             try
  //             {
  //               const auto fun =
  //                   sidesets_func
  //                       .cast<std::function<int(const RowVectorNd
  //                       &)>>();
  //               state.mesh->compute_boundary_ids(fun);
  //               return true;
  //             }
  //             catch (...)
  //             {
  //               {
  //               }
  //             }
  //             try
  //             {
  //               const auto fun = sidesets_func.cast<
  //                   std::function<int(const RowVectorNd &,
  //                   bool)>>();
  //               state.mesh->compute_boundary_ids(fun);
  //               return true;
  //             }
  //             catch (...)
  //             {
  //             }

  //             try
  //             {
  //               const auto fun = sidesets_func.cast<
  //                   std::function<int(const std::vector<int> &, bool)>>();
  //               state.mesh->compute_boundary_ids(fun);
  //               return true;
  //             }
  //             catch (...)
  //             {
  //             }

  //             throw pybind11::value_error(
  //                 "sidesets_func has invalid type, should be a function
  //                 (p)->int, (p, bool)->int, ([], bool)->int");
  //           }
  //         }();

  //         state.stats.compute_mesh_stats(*state.mesh);

  //         state.build_basis();

  //         state.assemble_rhs();
  //         state.assemble_mass_mat();
  //         state.solve_problem();

  //         return res;
  //       },
  //       "single solve function", py::kw_only(),
  //       py::arg("vertices") = Eigen::MatrixXd(),
  //       py::arg("cells") = Eigen::MatrixXi(),
  //       py::arg("sidesets_func") = py::none(), py::arg("log_level") = 2);
}
