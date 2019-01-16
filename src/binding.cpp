#include <polyfem/State.hpp>
#include <polyfem/Logger.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>

namespace py = pybind11;


int add(int i, int j) {
	return i + j;
}

PYBIND11_MODULE(polyfempy, m) {
	//TODO:
	// load mesh with lambda and V, F
	// out of vis mesh and evals of sol and pressure
	// load rhs form file or numpy

	py::class_<polyfem::State>(m, "Solver")
	.def(py::init<>())
	.def("load_parameters", [](polyfem::State &s, const std::string &json) {
		s.init(json::parse(json));
	}, 																"load PDE and problem parameters from a json")

	.def("set_log_level", [](polyfem::State &s, int log_level) {
		log_level = std::max(0, std::min(6, log_level));
		spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
	},																"sets polyfem log level, valid value between 0 (all logs) and 6 (no logs)")

	.def("load_mesh", [](polyfem::State &s) {
		s.load_mesh();
	},																"Loads a mesh from the 'mesh' field of the json and 'bc_tag' if any bc tags")
	.def("load_mesh", [](polyfem::State &s, std::string &path) {
		s.args["mesh"] = path;
		s.load_mesh();
	},																"Loads a mesh from the path and 'bc_tag' from the json if any bc tags")
	.def("load_mesh", [](polyfem::State &s, std::string &path, std::string &bc_tag) {
		s.args["mesh"] = path;
		s.args["bc_tag"] = bc_tag;
		s.load_mesh();
	},																"Loads a mesh and bc_tags from path")

	.def("set_rhs", [](polyfem::State &s, std::string &path) {
		s.args["rhs_path"] = path;
	}, 																"Loads the rhs from a file")
	.def("set_rhs", [](polyfem::State &s, const Eigen::MatrixXd &rhs) {
		s.rhs_in = rhs;
	}, 																"Sets the rhs")



	.def("solve", 			[](polyfem::State &s) {
		s.compute_mesh_stats();

		s.build_basis();
		s.build_polygonal_basis();

		s.assemble_rhs();
		s.assemble_stiffness_mat();

		s.solve_problem();
	}, 																"solve the pde")

	.def("compute_errors", 	&polyfem::State::compute_errors, 		"compute the error")

	.def("get_log",	[](polyfem::State &s) {
		std::stringstream ss;
		s.save_json(ss);
		return ss.str();
	},																"gets the log as json")
	.def("export_data", 	&polyfem::State::export_data,			"exports all data specified in the json")
	.def("export_vtu",	 	&polyfem::State::save_vtu,				"exports the solution as vtu")
	.def("export_wire", 	&polyfem::State::save_wire,				"exports wireframe of the mesh")


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
	}, 																					"returns the solution on a densly sampled mesh, use 'vismesh_rel_area' to control density",
	py::arg("boundary_only") = bool(false));

}