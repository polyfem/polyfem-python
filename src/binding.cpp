#include <polyfem/State.hpp>
#include <polyfem/AssemblerUtils.hpp>
#include <polyfem/Logger.hpp>
#include <polyfem/MeshUtils.hpp>

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>

namespace py = pybind11;

class ScalarAssemblers { };
class TensorAssemblers { };



namespace {
	void init_globals()
	{
		static bool initialized = false;

		if(!initialized)
		{
#ifndef WIN32
			setenv("GEO_NO_SIGNAL_HANDLER", "1", 1);
#endif

			GEO::initialize();

#ifdef USE_TBB
			const size_t MB = 1024*1024;
			const size_t stack_size = 64 * MB;
			unsigned int num_threads = std::max(1u, std::thread::hardware_concurrency());
			tbb::task_scheduler_init scheduler(num_threads, stack_size);
#endif

    		// Import standard command line arguments, and custom ones
			GEO::CmdLine::import_arg_group("standard");
			GEO::CmdLine::import_arg_group("pre");
			GEO::CmdLine::import_arg_group("algo");

			initialized = true;
		}
	}

}

PYBIND11_MODULE(polyfempy, m) {
	const auto &sa = py::class_<ScalarAssemblers>(m, "ScalarFormulations");
	for(auto &a : polyfem::AssemblerUtils::instance().scalar_assemblers())
		sa.attr(a.c_str()) = a;

	const auto &ta = py::class_<TensorAssemblers>(m, "TensorFormulations");
	for(auto &a : polyfem::AssemblerUtils::instance().tensor_assemblers())
		ta.attr(a.c_str()) = a;

	const auto &solver = py::class_<polyfem::State>(m, "Solver")

	.def(py::init<>())

	.def("settings", [](polyfem::State &self, const std::string &json) {
		init_globals();
		self.init(json::parse(json));
	},
	"load PDE and problem parameters from the settings",
	py::arg("json"))

	.def("set_log_level", [](polyfem::State &s, int log_level) {
		init_globals();
		log_level = std::max(0, std::min(6, log_level));
		spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
	},
	"sets polyfem log level, valid value between 0 (all logs) and 6 (no logs)",
	py::arg("log_level"))

	.def("load_mesh", [](polyfem::State &s) {
		init_globals();
		s.load_mesh();
	},
	"Loads a mesh from the 'mesh' field of the json and 'bc_tag' if any bc tags")

	.def("load_mesh", [](polyfem::State &s, const std::string &path) {
		init_globals();
		s.args["mesh"] = path;
		s.load_mesh();
	},
	"Loads a mesh from the path and 'bc_tag' from the json if any bc tags",
	py::arg("path"))
	.def("load_mesh", [](polyfem::State &s, const std::string &path, const std::string &bc_tag) {
		init_globals();
		s.args["mesh"] = path;
		s.args["bc_tag"] = bc_tag;
		s.load_mesh();
	},
	"Loads a mesh and bc_tags from path",
	py::arg("path"), py::arg("bc_tag_path"))
	.def("set_mesh", [](polyfem::State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F) {
		init_globals();

		GEO::Mesh M;
		if(V.cols() == 2)
			polyfem::to_geogram_mesh(V, F, M);
		else
			polyfem::to_geogram_mesh_3d(V, F, M);
		s.load_mesh(M, [](const polyfem::RowVectorNd&){ return -1; }, true);
	},
	"Loads a mesh from vertices and connectivity",
	py::arg("vertices"), py::arg("connectivity"))


	.def("set_boundary_side_set_from_bary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd&)> &boundary_marker) {
		init_globals();
		s.mesh->compute_boundary_ids(boundary_marker);
	},
	"Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face)",
	py::arg("boundary_marker"))
	.def("set_boundary_side_set_from_bary_and_boundary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd&, bool)> &boundary_marker) {
		init_globals();
		s.mesh->compute_boundary_ids(boundary_marker);
	},
	"Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face) and a flag that says if the element is boundary",
	py::arg("boundary_marker"))
	.def("set_boundary_side_set_from_v_ids", [](polyfem::State &s, const std::function<int(const std::vector<int>&, bool)> &boundary_marker) {
		init_globals();
		s.mesh->compute_boundary_ids(boundary_marker);
	},
	"Sets the side set for the boundary conditions, the functions takes the sorted list of vertex id and a flag that says if the element is boundary",
	py::arg("boundary_marker"))


	.def("set_rhs", [](polyfem::State &s, std::string &path) {
		init_globals();
		s.args["rhs_path"] = path;
	},
	"Loads the rhs from a file",
	py::arg("path"))
	.def("set_rhs", [](polyfem::State &s, const Eigen::MatrixXd &rhs) {
		init_globals();
		s.rhs_in = rhs;
	},
	"Sets the rhs",
	py::arg("matrix"))



	.def("solve",[](polyfem::State &s) {
		init_globals();

		s.compute_mesh_stats();

		s.build_basis();
		s.build_polygonal_basis();

		s.assemble_rhs();
		s.assemble_stiffness_mat();

		s.solve_problem();
	},
	"solve the pde")

	.def("compute_errors",
		&polyfem::State::compute_errors,
		"compute the error")

	.def("get_log",	[](polyfem::State &s) {
		std::stringstream ss;
		s.save_json(ss);
		return ss.str();
	},
	"gets the log as json")

	.def("export_data", 	&polyfem::State::export_data,			"exports all data specified in the settings")
	.def("export_vtu",	 	&polyfem::State::save_vtu,				"exports the solution as vtu", py::arg("path"))
	.def("export_wire", 	&polyfem::State::save_wire,				"exports wireframe of the mesh", py::arg("path"), py::arg("isolines") = false)


	.def("get_solution", 			[](const polyfem::State &s) { return s.sol;}, 		"returns the solution")
	.def("get_pressure", 			[](const polyfem::State &s) { return s.pressure;}, 	"returns the pressure")
	.def("get_sampled_solution", 	[](polyfem::State &s,  bool boundary_only) {
		Eigen::MatrixXd points;
		Eigen::MatrixXi tets;
		Eigen::MatrixXd discr;
		Eigen::MatrixXd fun;

		s.build_vis_mesh(points, tets, discr);


		s.interpolate_function(points.rows(), s.sol, fun, boundary_only);
		return  py::make_tuple(points, tets, fun);
	},
	"returns the solution on a densly sampled mesh, use 'vismesh_rel_area' to control density",
	py::arg("boundary_only") = bool(false))

	.def("get_stresses", 	[](polyfem::State &s,  bool boundary_only) {
		Eigen::MatrixXd points;
		Eigen::MatrixXi tets;
		Eigen::MatrixXd discr;
		Eigen::MatrixXd fun;

		s.build_vis_mesh(points, tets, discr);
		s.compute_tensor_value(points.rows(), s.sol, fun, boundary_only);

		return fun;
	},
	"returns the stress tensor on a densly sampled mesh, use 'vismesh_rel_area' to control density",
	py::arg("boundary_only") = bool(false))

	.def("get_sampled_mises", 	[](polyfem::State &s,  bool boundary_only) {
		Eigen::MatrixXd points;
		Eigen::MatrixXi tets;
		Eigen::MatrixXd discr;
		Eigen::MatrixXd fun;

		s.build_vis_mesh(points, tets, discr);
		s.compute_scalar_value(points.rows(), s.sol, fun, boundary_only);

		return fun;
	},
	"returns the von mises stresses on a densly sampled mesh, use 'vismesh_rel_area' to control density",
	py::arg("boundary_only") = bool(false))

	.def("get_sampled_mises_avg", 	[](polyfem::State &s,  bool boundary_only) {
		Eigen::MatrixXd points;
		Eigen::MatrixXi tets;
		Eigen::MatrixXd discr;
		Eigen::MatrixXd fun, tfun;

		s.build_vis_mesh(points, tets, discr);
		s.average_grad_based_function(points.rows(), s.sol, fun, tfun, boundary_only);

		return py::make_tuple(fun, tfun);
	},
	"returns the von mises stresses and stress tensor averaged around a vertex on a densly sampled mesh, use 'vismesh_rel_area' to control density",
	py::arg("boundary_only") = bool(false));

	solver.doc() = "Polyfem solver";

}