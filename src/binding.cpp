#include <polyfem/State.hpp>
#include <polyfem/AssemblerUtils.hpp>
#include <polyfem/Logger.hpp>
#include <polyfem/MeshUtils.hpp>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <igl/remove_duplicate_vertices.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

	namespace py = pybind11;

class ScalarAssemblers
{
};
class TensorAssemblers
{
};

class PDEs
{ };

//TODO add save_time_sequence

namespace
{
void init_globals(polyfem::State &state)
{
	static bool initialized = false;

	if (!initialized)
	{
#ifndef WIN32
		setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

		GEO::initialize();

#ifdef USE_TBB
		const size_t MB = 1024 * 1024;
		const size_t stack_size = 64 * MB;
		unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
		tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

		// Import standard command line arguments, and custom ones
		GEO::CmdLine::import_arg_group("standard");
		GEO::CmdLine::import_arg_group("pre");
		GEO::CmdLine::import_arg_group("algo");

		state.init_logger(std::cout, 2);

		initialized = true;
	}
}

} // namespace

PYBIND11_MODULE(polyfempy, m)
{
	const auto &pdes = py::class_<PDEs>(m, "PDEs");

	const auto &sa = py::class_<ScalarAssemblers>(m, "ScalarFormulations");
	for (auto &a : polyfem::AssemblerUtils::instance().scalar_assemblers()){
		sa.attr(a.c_str()) = a;
		pdes.attr(a.c_str()) = a;
	}

	const auto &ta = py::class_<TensorAssemblers>(m, "TensorFormulations");
	for (auto &a : polyfem::AssemblerUtils::instance().tensor_assemblers()){
		ta.attr(a.c_str()) = a;
		pdes.attr(a.c_str()) = a;
	}

	ta.attr("NonLinearElasticity") = "NonLinearElasticity";
	pdes.attr("NonLinearElasticity") = "NonLinearElasticity";

	pdes.doc() = "List of supported partial differential equations";

	m.def("is_tensor", [](const std::string &pde) {
		const auto &assembler = polyfem::AssemblerUtils::instance();
		return assembler.is_tensor(pde);
	},
			"returns true if the pde is tensorial", py::arg("pde"));

	const auto setting_lambda = [](polyfem::State &self, const py::object &json) {
		init_globals(self);
		py::scoped_ostream_redirect output;
		const std::string json_string = py::str(json);
		self.init(json::parse(json_string));
	};
	auto &solver = py::class_<polyfem::State>(m, "Solver")
					   .def(py::init<>())

					   .def("settings", setting_lambda,
							"load PDE and problem parameters from the settings", py::arg("json"))

					   .def("set_settings", setting_lambda,
							"load PDE and problem parameters from the settings", py::arg("json"))

					   .def("set_log_level", [](polyfem::State &s, int log_level) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   log_level = std::max(0, std::min(6, log_level));
						   spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
					   },
							"sets polyfem log level, valid value between 0 (all logs) and 6 (no logs)", py::arg("log_level"))

					   .def("load_mesh_from_settings", [](polyfem::State &s, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.args["normalize_mesh"] = normalize_mesh;
						   s.args["n_refs"] = n_refs;
						   s.args["boundary_id_threshold"] = boundary_id_threshold;
						   s.args["vismesh_rel_area"] = vismesh_rel_area;
						   s.load_mesh();
					   },
							"Loads a mesh from the 'mesh' field of the json and 'bc_tag' if any bc tags", py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def("load_mesh_from_path", [](polyfem::State &s, const std::string &path, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.args["mesh"] = path;
						   s.args["normalize_mesh"] = normalize_mesh;
						   s.args["n_refs"] = n_refs;
						   s.args["boundary_id_threshold"] = boundary_id_threshold;
						   s.args["vismesh_rel_area"] = vismesh_rel_area;
						   s.load_mesh();
					   },
							"Loads a mesh from the path and 'bc_tag' from the json if any bc tags", py::arg("path"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def("load_mesh_from_path_and_tags", [](polyfem::State &s, const std::string &path, const std::string &bc_tag, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.args["mesh"] = path;
						   s.args["bc_tag"] = bc_tag;
						   s.args["normalize_mesh"] = normalize_mesh;
						   s.args["n_refs"] = n_refs;
						   s.args["boundary_id_threshold"] = boundary_id_threshold;
						   s.args["vismesh_rel_area"] = vismesh_rel_area;
						   s.load_mesh();
					   },
							"Loads a mesh and bc_tags from path", py::arg("path"), py::arg("bc_tag_path"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def("set_mesh", [](polyfem::State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;

						   if (V.cols() == 2)
							   s.mesh = std::make_unique<polyfem::Mesh2D>();
						   else
							   s.mesh = std::make_unique<polyfem::Mesh3D>();
						   s.mesh->build_from_matrices(V, F);

						   s.args["normalize_mesh"] = normalize_mesh;
						   s.args["n_refs"] = n_refs;
						   s.args["boundary_id_threshold"] = boundary_id_threshold;
						   s.args["vismesh_rel_area"] = vismesh_rel_area;

						   s.load_mesh();
					   },
							"Loads a mesh from vertices and connectivity", py::arg("vertices"), py::arg("connectivity"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def("set_high_order_mesh", [](polyfem::State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &nodes_pos, const std::vector<std::vector<int>> &nodes_indices, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;

						   if (V.cols() == 2)
							   s.mesh = std::make_unique<polyfem::Mesh2D>();
						   else
							   s.mesh = std::make_unique<polyfem::Mesh3D>();
						   s.mesh->build_from_matrices(V, F);
						   s.mesh->attach_higher_order_nodes(nodes_pos, nodes_indices);

						   s.args["normalize_mesh"] = normalize_mesh;
						   s.args["n_refs"] = n_refs;
						   s.args["boundary_id_threshold"] = boundary_id_threshold;
						   s.args["vismesh_rel_area"] = vismesh_rel_area;

						   s.load_mesh();
					   },
							"Loads an high order mesh from vertices, connectivity, nodes, and node indices mapping element to nodes", py::arg("vertices"), py::arg("connectivity"), py::arg("nodes_pos"), py::arg("nodes_indices"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def("set_boundary_side_set_from_bary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd &)> &boundary_marker) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.mesh->compute_boundary_ids(boundary_marker);
					   },
							"Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face)", py::arg("boundary_marker"))
					   .def("set_boundary_side_set_from_bary_and_boundary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd &, bool)> &boundary_marker) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.mesh->compute_boundary_ids(boundary_marker);
					   },
							"Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face) and a flag that says if the element is boundary", py::arg("boundary_marker"))
					   .def("set_boundary_side_set_from_v_ids", [](polyfem::State &s, const std::function<int(const std::vector<int> &, bool)> &boundary_marker) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.mesh->compute_boundary_ids(boundary_marker);
					   },
							"Sets the side set for the boundary conditions, the functions takes the sorted list of vertex id and a flag that says if the element is boundary", py::arg("boundary_marker"))

					   .def("set_rhs_from_path", [](polyfem::State &s, std::string &path) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.args["rhs_path"] = path;
					   },
							"Loads the rhs from a file", py::arg("path"))
					   .def("set_rhs", [](polyfem::State &s, const Eigen::MatrixXd &rhs) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.rhs_in = rhs;
					   },
							"Sets the rhs", py::arg("matrix"))

					   .def("solve", [](polyfem::State &s) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;
						   s.compute_mesh_stats();

						   s.build_basis();

						   s.assemble_rhs();
						   s.assemble_stiffness_mat();

						   s.solve_export_to_file = false;
						   s.solution_frames.clear();
						   s.solve_problem();
						   s.solve_export_to_file = true;
					   },
							"solve the pde")

					   .def("compute_errors", [](polyfem::State &s) {
						   init_globals(s);
						   py::scoped_ostream_redirect output;

						   s.compute_errors();
					   },
							"compute the error")

					   .def("get_log", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   std::stringstream ss;
						   s.save_json(ss);
						   return ss.str();
					   },
							"gets the log as json")

					   .def("export_data", [](polyfem::State &s) { py::scoped_ostream_redirect output; s.export_data(); }, "exports all data specified in the settings")
					   .def("export_vtu", [](polyfem::State &s, std::string &path, bool boundary_only) {
						   py::scoped_ostream_redirect output;
						   const bool tmp = s.args["export"]["vis_boundary_only"];
						   s.args["export"]["vis_boundary_only"] = boundary_only;
						   s.save_vtu(path, 0);
						   s.args["export"]["vis_boundary_only"] = tmp;
						},
						"exports the solution as vtu", py::arg("path"), py::arg("boundary_only") = bool(false))
					   .def("export_wire", [](polyfem::State &s, std::string &path, bool isolines) { py::scoped_ostream_redirect output; s.save_wire(path, isolines); }, "exports wireframe of the mesh", py::arg("path"), py::arg("isolines") = false)

					   .def("get_solution", [](const polyfem::State &s) { return s.sol; }, "returns the solution")
					   .def("get_pressure", [](const polyfem::State &s) { return s.pressure; }, "returns the pressure")
					   .def("get_sampled_solution", [](polyfem::State &s, bool boundary_only) {
						   py::scoped_ostream_redirect output;
						   Eigen::MatrixXd points;
						   Eigen::MatrixXi tets;
						   Eigen::MatrixXd discr;
						   Eigen::MatrixXd fun;

						   const bool tmp = s.args["export"]["vis_boundary_only"];
						   s.args["export"]["vis_boundary_only"] = boundary_only;

						   s.build_vis_mesh(points, tets, discr);

						   s.interpolate_function(points.rows(), s.sol, fun, boundary_only);

						   s.args["export"]["vis_boundary_only"] = tmp;

						   return py::make_tuple(points, tets, fun);
					   },
							"returns the solution on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def("get_stresses", [](polyfem::State &s, bool boundary_only) {
						   py::scoped_ostream_redirect output;
						   Eigen::MatrixXd points;
						   Eigen::MatrixXi tets;
						   Eigen::MatrixXd discr;
						   Eigen::MatrixXd fun;

						   const bool tmp = s.args["export"]["vis_boundary_only"];
						   s.args["export"]["vis_boundary_only"] = boundary_only;

						   s.build_vis_mesh(points, tets, discr);
						   s.compute_tensor_value(points.rows(), s.sol, fun, boundary_only);

						   s.args["export"]["vis_boundary_only"] = tmp;

						   return fun;
					   },
							"returns the stress tensor on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def("get_sampled_mises", [](polyfem::State &s, bool boundary_only) {
						   py::scoped_ostream_redirect output;
						   Eigen::MatrixXd points;
						   Eigen::MatrixXi tets;
						   Eigen::MatrixXd discr;
						   Eigen::MatrixXd fun;

						   const bool tmp = s.args["export"]["vis_boundary_only"];
						   s.args["export"]["vis_boundary_only"] = boundary_only;

						   s.build_vis_mesh(points, tets, discr);
						   s.compute_scalar_value(points.rows(), s.sol, fun, boundary_only);

						   s.args["export"]["vis_boundary_only"] = tmp;

						   return fun;
					   },
							"returns the von mises stresses on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def("get_sampled_mises_avg", [](polyfem::State &s, bool boundary_only) {
						   py::scoped_ostream_redirect output;
						   Eigen::MatrixXd points;
						   Eigen::MatrixXi tets;
						   Eigen::MatrixXd discr;
						   Eigen::MatrixXd fun, tfun;

						   const bool tmp = s.args["export"]["vis_boundary_only"];
						   s.args["export"]["vis_boundary_only"] = boundary_only;

						   s.build_vis_mesh(points, tets, discr);
						   s.average_grad_based_function(points.rows(), s.sol, fun, tfun, boundary_only);

						   s.args["export"]["vis_boundary_only"] = tmp;

						   return py::make_tuple(fun, tfun);
					   },
							"returns the von mises stresses and stress tensor averaged around a vertex on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def("get_sampled_traction_forces", [](polyfem::State &s, const bool apply_displacement, const bool compute_avg) {
						   py::scoped_ostream_redirect output;

						   if (!s.mesh)
							   throw pybind11::value_error("Load the mesh first!");
						   if (!s.mesh->is_volume())
							   throw pybind11::value_error("This function works only on volumetric meshes!");
						   if (s.problem->is_scalar())
							   throw pybind11::value_error("Define a tensor problem!");




						   Eigen::MatrixXd result;

						   Eigen::MatrixXd v_surf;
						   Eigen::MatrixXi f_surf;
						   const double epsilon = 1e-10;

						   {
							   const polyfem::Mesh3D &mesh3d = *dynamic_cast<polyfem::Mesh3D *>(s.mesh.get());
							   Eigen::MatrixXd points(mesh3d.n_vertices(), 3);
							   Eigen::MatrixXi tets(mesh3d.n_cells(), 4);

							   for(int t = 0; t < mesh3d.n_cells(); ++t)
							   {
								   if (mesh3d.n_cell_vertices(t) != 4)
									   throw pybind11::value_error("Works only with tet meshes!");

									for(int i = 0; i < 4; ++i)
								   		tets(t, i) = mesh3d.cell_vertex(t, i);
							   }

							   for(int p = 0; p < mesh3d.n_vertices(); ++p)
								   points.row(p) << mesh3d.point(p);

							   Eigen::MatrixXi f_surf_tmp, _;
							   igl::boundary_facets(tets, f_surf_tmp);
							   igl::remove_unreferenced(points, f_surf_tmp, v_surf, f_surf, _);
						   }

						   if(apply_displacement)
							   s.interpolate_boundary_tensor_function(v_surf, f_surf, s.sol, s.sol, compute_avg, result);
						   else
							   s.interpolate_boundary_tensor_function(v_surf, f_surf, s.sol, compute_avg, result);

						   return py::make_tuple(v_surf, f_surf, result);
					   },
							"returns the traction forces computed on the surface", py::arg("apply_displacement") = bool(false), py::arg("compute_avg") = bool(true))

					   ////////////////////////////////////////////////////////////////////////////////////////////
					   .def("get_sampled_points_frames", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   assert(!s.solution_frames.empty());

						   std::vector<Eigen::MatrixXd> pts;

						   for (const auto &sol : s.solution_frames)
						   {
							   pts.push_back(sol.points);
						   }

						   return pts;
					   },
							"returns the points frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def("get_sampled_connectivity_frames", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   assert(!s.solution_frames.empty());

						   std::vector<Eigen::MatrixXi> tets;

						   for (const auto &sol : s.solution_frames)
							   tets.push_back(sol.connectivity);

						   return tets;
					   },
							"returns the connectivity frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def("get_sampled_solution_frames", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   assert(!s.solution_frames.empty());

						   std::vector<Eigen::MatrixXd> fun;

						   for (const auto &sol : s.solution_frames)
						   {
							   fun.push_back(sol.solution);
						   }

						   return fun;
					   },
							"returns the solution frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def("get_sampled_mises_frames", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   assert(!s.solution_frames.empty());

						   std::vector<Eigen::MatrixXd> mises;

						   for (const auto &sol : s.solution_frames)
							   mises.push_back(sol.scalar_value);

						   return mises;
					   },
							"returns the von mises stresses frames on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def("get_sampled_mises_avg_frames", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   assert(!s.solution_frames.empty());

						   std::vector<Eigen::MatrixXd> mises;

						   for (const auto &sol : s.solution_frames)
							   mises.push_back(sol.scalar_value_avg);

						   return mises;
					   },
							"returns the von mises stresses per frame averaged around a vertex on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def("get_boundary_sidesets", [](polyfem::State &s) {
						   py::scoped_ostream_redirect output;
						   Eigen::MatrixXd points;
						   Eigen::MatrixXi faces;
						   Eigen::MatrixXd sidesets;

						   s.get_sidesets(points, faces, sidesets);

						   return py::make_tuple(points, faces, sidesets);
					   },
							"exports get the boundary sideset, edges in 2d or trangles in 3d");

	solver.doc() = "Polyfem solver";
}