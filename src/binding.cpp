#include <polyfem/State.hpp>
#include <polyfem/AssemblerUtils.hpp>
#include <polyfem/Logger.hpp>
#include <polyfem/MeshUtils.hpp>
#include <polyfem/GenericProblem.hpp>
#include <polyfem/StringUtils.hpp>

#include "raster.hpp"

#include <geogram/basic/command_line.h>
#include <geogram/basic/command_line_args.h>

#include <igl/remove_duplicate_vertices.h>
#include <igl/boundary_facets.h>
#include <igl/remove_unreferenced.h>

#ifdef USE_TBB
#include <tbb/task_scheduler_init.h>
#include <thread>
#endif

#include <pybind11_json/pybind11_json.hpp>

#include <pybind11/pybind11.h>
#include <pybind11/eigen.h>
#include <pybind11/functional.h>
#include <pybind11/stl.h>
#include <pybind11/iostream.h>

namespace py = pybind11;

typedef std::function<Eigen::MatrixXd(double x, double y, double z)> BCFuncV;
typedef std::function<double(double x, double y, double z)> BCFuncS;

class ScalarAssemblers
{
};
class TensorAssemblers
{
};

class PDEs
{
};

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
	for (auto &a : polyfem::AssemblerUtils::scalar_assemblers())
	{
		sa.attr(a.c_str()) = a;
		pdes.attr(a.c_str()) = a;
	}

	const auto &ta = py::class_<TensorAssemblers>(m, "TensorFormulations");
	for (auto &a : polyfem::AssemblerUtils::tensor_assemblers())
	{
		ta.attr(a.c_str()) = a;
		pdes.attr(a.c_str()) = a;
	}

	ta.attr("NonLinearElasticity") = "NonLinearElasticity";
	pdes.attr("NonLinearElasticity") = "NonLinearElasticity";

	pdes.doc() = "List of supported partial differential equations";

	m.def(
		"is_tensor", [](const std::string &pde)
		{ return polyfem::AssemblerUtils::is_tensor(pde); },
		"returns true if the pde is tensorial", py::arg("pde"));

	const auto setting_lambda = [](polyfem::State &self, const py::object &settings)
	{
		using namespace polyfem;

		init_globals(self);
		// py::scoped_ostream_redirect output;
		const std::string json_string = py::str(settings);
		self.init(json::parse(json_string));
	};

	const auto rendering_lambda = [] (const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &faces, int width, int height, const Eigen::MatrixXd &camera_positionm, const double camera_fov, const double camera_near, const double camera_far, const bool is_perspective, const Eigen::MatrixXd &lookatm, const Eigen::MatrixXd &upm, const Eigen::MatrixXd &ambient_lightm, const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> &lights, const Eigen::MatrixXd &diffuse_colorm, const Eigen::MatrixXd &specular_colorm, const double specular_exponent)
		{
			using namespace renderer;
			using namespace Eigen;

			Eigen::Vector3d camera_position(0, 0, 1);
			Eigen::Vector3d lookat(0, 0, 0);
			Eigen::Vector3d up(0, 0, 1);
			Eigen::Vector3d ambient_light(0.1, 0.1, 0.1);
			Eigen::Vector3d diffuse_color(1, 0, 0);
			Eigen::Vector3d specular_color(1, 0, 0);

			if (camera_positionm.size() > 0 && camera_positionm.size() != 3)
				throw pybind11::value_error("camera_position have size 3");
			if (camera_positionm.size() == 3)
				camera_position << camera_positionm(0), camera_positionm(1), camera_positionm(2);
			if (lookatm.size() > 0 && lookatm.size() != 3)
				throw pybind11::value_error("lookat have size 3");
			if (lookatm.size() == 3)
				lookat << lookatm(0), lookatm(1), lookatm(2);
			if (upm.size() > 0 && upm.size() != 3)
				throw pybind11::value_error("up have size 3");
			if (upm.size() == 3)
				up << upm(0), upm(1), upm(2);
			if (ambient_lightm.size() > 0 && ambient_lightm.size() != 3)
				throw pybind11::value_error("ambient_light have size 3");
			if (ambient_lightm.size() == 3)
				ambient_light << ambient_lightm(0), ambient_lightm(1), ambient_lightm(2);
			if (diffuse_colorm.size() > 0 && diffuse_colorm.size() != 3)
				throw pybind11::value_error("diffuse_color have size 3");
			if (diffuse_colorm.size() == 3)
				diffuse_color << diffuse_colorm(0), diffuse_colorm(1), diffuse_colorm(2);
			if (specular_colorm.size() > 0 && specular_colorm.size() != 3)
				throw pybind11::value_error("specular_color have size 3");
			if (specular_colorm.size() == 3)
				specular_color << specular_colorm(0), specular_colorm(1), specular_colorm(2);

			Material material;
			material.diffuse_color = diffuse_color;
			material.specular_color = specular_color;
			material.specular_exponent = specular_exponent;

			Eigen::Matrix<FrameBufferAttributes, Eigen::Dynamic, Eigen::Dynamic> frameBuffer(width, height);
			UniformAttributes uniform;

			const Vector3d gaze = lookat - camera_position;
			const Vector3d w = -gaze.normalized();
			const Vector3d u = up.cross(w).normalized();
			const Vector3d v = w.cross(u);

			Matrix4d M_cam_inv;
			M_cam_inv << u(0), v(0), w(0), camera_position(0),
				u(1), v(1), w(1), camera_position(1),
				u(2), v(2), w(2), camera_position(2),
				0, 0, 0, 1;

			uniform.M_cam = M_cam_inv.inverse();

			{
				const double camera_ar = double(width) / height;
				const double tan_angle = tan(camera_fov / 2);
				const double n = - camera_near;
				const double f = - camera_far;
				const double t = tan_angle * n;
				const double b = -t;
				const double r = t * camera_ar;
				const double l = -r;

				uniform.M_orth << 2 / (r - l), 0, 0, -(r + l) / (r - l),
					0, 2 / (t - b), 0, -(t + b) / (t - b),
					0, 0, 2 / (n - f), -(n + f) / (n - f),
					0, 0, 0, 1;
				Matrix4d P;
				if (is_perspective)
				{
					P << n, 0, 0, 0,
						0, n, 0, 0,
						0, 0, n + f, -f * n,
						0, 0, 1, 0;
				}
				else
				{
					P << 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1;
				}

				uniform.M = uniform.M_orth * P * uniform.M_cam;
			}

			Program program;
			program.VertexShader = [&](const VertexAttributes &va, const UniformAttributes &uniform)
			{
				VertexAttributes out;
				out.position = uniform.M * va.position;
				Vector3d color = ambient_light;

				Vector3d hit(va.position(0), va.position(1), va.position(2));
				for (const auto &l : lights)
				{
					Vector3d Li = (l.first - hit).normalized();
					Vector3d N = va.normal;
					Vector3d diffuse = va.material.diffuse_color * std::max(Li.dot(N), 0.0);
					Vector3d H;
					if (is_perspective)
					{
						H = (Li + (camera_position - hit).normalized()).normalized();
					}
					else
					{
						H = (Li - gaze.normalized()).normalized();
					}
					const Vector3d specular = va.material.specular_color * std::pow(std::max(N.dot(H), 0.0), va.material.specular_exponent);
					const Vector3d D = l.first - hit;
					color += (diffuse + specular).cwiseProduct(l.second) / D.squaredNorm();
				}
				out.color = color;
				return out;
			};

			program.FragmentShader = [](const VertexAttributes &va, const UniformAttributes &uniform)
			{
				FragmentAttributes out(va.color(0), va.color(1), va.color(2));
				out.depth = -va.position(2);
				return out;
			};

			program.BlendingShader = [](const FragmentAttributes &fa, const FrameBufferAttributes &previous)
			{
				if (fa.depth < previous.depth)
				{
					FrameBufferAttributes out(fa.color[0] * 255, fa.color[1] * 255, fa.color[2] * 255, fa.color[3] * 255);
					out.depth = fa.depth;
					return out;
				}
				else
				{
					return previous;
				}
			};

			Eigen::MatrixXd vnormals(vertices.rows(), 3);
			//    Eigen::MatrixXd areas(tmp.rows(), 1);
			vnormals.setZero();
			//    areas.setZero();
			Eigen::MatrixXd fnormals(faces.rows(), 3);

			for (int i = 0; i < faces.rows(); ++i)
			{
				const Vector3d l1 = vertices.row(faces(i, 1)) - vertices.row(faces(i, 0));
				const Vector3d l2 = vertices.row(faces(i, 2)) - vertices.row(faces(i, 0));
				const auto nn = l1.cross(l2);
				const double area = nn.norm();
				fnormals.row(i) = nn / area;

				for (int j = 0; j < 3; j++)
				{
					int vid = faces(i, j);
					vnormals.row(vid) += nn;
					//    areas(vid) += area;
				}
			}

			std::vector<VertexAttributes> vertex_attributes;
			for (int i = 0; i < faces.rows(); ++i)
			{
				for (int j = 0; j < 3; j++)
				{
					int vid = faces(i, j);
					VertexAttributes va(vertices(vid, 0), vertices(vid, 1), vertices(vid, 2));
					va.material = material;
					va.normal = vnormals.row(vid).normalized();
					vertex_attributes.push_back(va);
				}
			}

			rasterize_triangles(program, uniform, vertex_attributes, frameBuffer);

			std::vector<uint8_t> image;
			framebuffer_to_uint8(frameBuffer, image);

			return image;
	};

	auto &solver = py::class_<polyfem::State>(m, "Solver")
					   .def(py::init<>())

					   .def("settings", setting_lambda,
							"load PDE and problem parameters from the settings", py::arg("json"))

					   .def("set_settings", setting_lambda,
							"load PDE and problem parameters from the settings", py::arg("json"))

					   .def(
						   "set_log_level", [](polyfem::State &s, int log_level)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   log_level = std::max(0, std::min(6, log_level));
							   spdlog::set_level(static_cast<spdlog::level::level_enum>(log_level));
						   },
						   "sets polyfem log level, valid value between 0 (all logs) and 6 (no logs)", py::arg("log_level"))

					   .def(
						   "load_mesh_from_settings", [](polyfem::State &s, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.args["normalize_mesh"] = normalize_mesh;
							   s.args["n_refs"] = n_refs;
							   s.args["boundary_id_threshold"] = boundary_id_threshold;
							   s.args["vismesh_rel_area"] = vismesh_rel_area;
							   s.load_mesh();
						   },
						   "Loads a mesh from the 'mesh' field of the json and 'bc_tag' if any bc tags", py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def(
						   "load_mesh_from_path", [](polyfem::State &s, const std::string &path, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.args["mesh"] = path;
							   s.args["normalize_mesh"] = normalize_mesh;
							   s.args["n_refs"] = n_refs;
							   s.args["boundary_id_threshold"] = boundary_id_threshold;
							   s.args["vismesh_rel_area"] = vismesh_rel_area;
							   s.load_mesh();
						   },
						   "Loads a mesh from the path and 'bc_tag' from the json if any bc tags", py::arg("path"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def(
						   "load_mesh_from_path_and_tags", [](polyfem::State &s, const std::string &path, const std::string &bc_tag, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.args["mesh"] = path;
							   s.args["bc_tag"] = bc_tag;
							   s.args["normalize_mesh"] = normalize_mesh;
							   s.args["n_refs"] = n_refs;
							   s.args["boundary_id_threshold"] = boundary_id_threshold;
							   s.args["vismesh_rel_area"] = vismesh_rel_area;
							   s.load_mesh();
						   },
						   "Loads a mesh and bc_tags from path", py::arg("path"), py::arg("bc_tag_path"), py::arg("normalize_mesh") = bool(false), py::arg("vismesh_rel_area") = double(0.00001), py::arg("n_refs") = int(0), py::arg("boundary_id_threshold") = double(-1))

					   .def(
						   "set_mesh", [](polyfem::State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;

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

					   .def(
						   "set_high_order_mesh", [](polyfem::State &s, const Eigen::MatrixXd &V, const Eigen::MatrixXi &F, const Eigen::MatrixXd &nodes_pos, const std::vector<std::vector<int>> &nodes_indices, const bool normalize_mesh, const double vismesh_rel_area, const int n_refs, const double boundary_id_threshold)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;

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

					   .def(
						   "set_boundary_side_set_from_bary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd &)> &boundary_marker)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.mesh->compute_boundary_ids(boundary_marker);
						   },
						   "Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face)", py::arg("boundary_marker"))
					   .def(
						   "set_boundary_side_set_from_bary_and_boundary", [](polyfem::State &s, const std::function<int(const polyfem::RowVectorNd &, bool)> &boundary_marker)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.mesh->compute_boundary_ids(boundary_marker);
						   },
						   "Sets the side set for the boundary conditions, the functions takes the barycenter of the boundary (edge or face) and a flag that says if the element is boundary", py::arg("boundary_marker"))
					   .def(
						   "set_boundary_side_set_from_v_ids", [](polyfem::State &s, const std::function<int(const std::vector<int> &, bool)> &boundary_marker)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.mesh->compute_boundary_ids(boundary_marker);
						   },
						   "Sets the side set for the boundary conditions, the functions takes the sorted list of vertex id and a flag that says if the element is boundary", py::arg("boundary_marker"))

					   .def(
						   "set_rhs_from_path", [](polyfem::State &s, std::string &path)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.args["rhs_path"] = path;
						   },
						   "Loads the rhs from a file", py::arg("path"))
					   .def(
						   "set_rhs", [](polyfem::State &s, const Eigen::MatrixXd &rhs)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.rhs_in = rhs;
						   },
						   "Sets the rhs", py::arg("matrix"))

					   .def(
						   "solve", [](polyfem::State &s)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
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
					   .def(
						   "init_timestepping", [](polyfem::State &s, const double t0, const double dt)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;
							   s.compute_mesh_stats();

							   s.build_basis();

							   s.assemble_rhs();
							   s.assemble_stiffness_mat();

							   s.solution_frames.clear();
							   Eigen::VectorXd _;
							   s.init_transient(_);
							   const polyfem::RhsAssembler &rhs_assembler = *s.step_data.rhs_assembler;
							   s.solve_transient_tensor_non_linear_init(t0, dt, rhs_assembler);
						   },
						   "initialize timestepping", py::arg("t0"), py::arg("dt"))
					   .def(
						   "step_in_time", [](polyfem::State &s, const double t0, const double dt, const int t)
						   {
							   //    std::cout << s.step_data.rhs_assembler->problem_ << std::endl;
							   //    std::cout << &s.step_data.nl_problem->rhs_assembler.problem_ << std::endl;

							   //    std::cout << &s.step_data.nl_problem->rhs_assembler << std::endl;
							   //    std::cout << s.step_data.rhs_assembler << std::endl;

							   //    std::cout << s.step_data.rhs_assembler->problem_->is_rhs_zero() << std::endl;
							   //    std::cout << s.step_data.nl_problem->rhs_assembler.problem_->is_rhs_zero() << std::endl;
							   json solver_info;
							   s.solve_transient_tensor_non_linear_step(t0, dt, t, solver_info);
							   return solver_info;
						   },
						   "step in time", py::arg("t0"), py::arg("dt"), py::arg("t"))

					   .def(
						   "compute_errors", [](polyfem::State &s)
						   {
							   init_globals(s);
							   //    py::scoped_ostream_redirect output;

							   s.compute_errors();
						   },
						   "compute the error")

					   .def(
						   "get_log", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   json out;
							   s.save_json(out);
							   return out;
						   },
						   "gets the log")

					   //    .def("export_data", [](polyfem::State &s) { py::scoped_ostream_redirect output; s.export_data(); }, "exports all data specified in the settings")
					   .def(
						   "export_vtu", [](polyfem::State &s, std::string &path, bool boundary_only)
						   {
							   //    py::scoped_ostream_redirect output;
							   const bool tmp = s.args["export"]["vis_boundary_only"];
							   s.args["export"]["vis_boundary_only"] = boundary_only;
							   s.save_vtu(path, 0);
							   s.args["export"]["vis_boundary_only"] = tmp;
						   },
						   "exports the solution as vtu", py::arg("path"), py::arg("boundary_only") = bool(false))
					   //    .def("export_wire", [](polyfem::State &s, std::string &path, bool isolines) { py::scoped_ostream_redirect output; s.save_wire(path, isolines); }, "exports wireframe of the mesh", py::arg("path"), py::arg("isolines") = false)

					   .def(
						   "get_solution", [](const polyfem::State &s)
						   { return s.sol; },
						   "returns the solution")
					   .def(
						   "get_pressure", [](const polyfem::State &s)
						   { return s.pressure; },
						   "returns the pressure")
					   .def(
						   "get_sampled_solution", [](polyfem::State &s, bool boundary_only)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi tets;
							   Eigen::MatrixXi el_id;
							   Eigen::MatrixXd discr;
							   Eigen::MatrixXd fun;

							   const bool tmp = s.args["export"]["vis_boundary_only"];
							   s.args["export"]["vis_boundary_only"] = boundary_only;

							   s.build_vis_mesh(points, tets, el_id, discr);

							   Eigen::MatrixXd ids(points.rows(), 1);
							   for (int i = 0; i < points.rows(); ++i)
							   {
								   ids(i) = s.mesh->get_body_id(el_id(i));
							   }

							   s.interpolate_function(points.rows(), s.sol, fun, boundary_only);

							   s.args["export"]["vis_boundary_only"] = tmp;

							   return py::make_tuple(points, tets, el_id, ids, fun);
						   },
						   "returns the solution on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def(
						   "get_stresses", [](polyfem::State &s, bool boundary_only)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi tets;
							   Eigen::MatrixXi el_id;
							   Eigen::MatrixXd discr;
							   Eigen::MatrixXd fun;

							   const bool tmp = s.args["export"]["vis_boundary_only"];
							   s.args["export"]["vis_boundary_only"] = boundary_only;

							   s.build_vis_mesh(points, tets, el_id, discr);
							   s.compute_tensor_value(points.rows(), s.sol, fun, boundary_only);

							   s.args["export"]["vis_boundary_only"] = tmp;

							   return fun;
						   },
						   "returns the stress tensor on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def(
						   "get_sampled_mises", [](polyfem::State &s, bool boundary_only)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi tets;
							   Eigen::MatrixXi el_id;
							   Eigen::MatrixXd discr;
							   Eigen::MatrixXd fun;

							   const bool tmp = s.args["export"]["vis_boundary_only"];
							   s.args["export"]["vis_boundary_only"] = boundary_only;

							   s.build_vis_mesh(points, tets, el_id, discr);
							   s.compute_scalar_value(points.rows(), s.sol, fun, boundary_only);

							   s.args["export"]["vis_boundary_only"] = tmp;

							   return fun;
						   },
						   "returns the von mises stresses on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def(
						   "get_sampled_mises_avg", [](polyfem::State &s, bool boundary_only)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi tets;
							   Eigen::MatrixXi el_id;
							   Eigen::MatrixXd discr;
							   Eigen::MatrixXd fun, tfun;

							   const bool tmp = s.args["export"]["vis_boundary_only"];
							   s.args["export"]["vis_boundary_only"] = boundary_only;

							   s.build_vis_mesh(points, tets, el_id, discr);
							   s.average_grad_based_function(points.rows(), s.sol, fun, tfun, boundary_only);

							   s.args["export"]["vis_boundary_only"] = tmp;

							   return py::make_tuple(fun, tfun);
						   },
						   "returns the von mises stresses and stress tensor averaged around a vertex on a densly sampled mesh, use 'vismesh_rel_area' to control density", py::arg("boundary_only") = bool(false))

					   .def(
						   "get_sampled_traction_forces", [](polyfem::State &s, const bool apply_displacement, const bool compute_avg)
						   {
							   //    py::scoped_ostream_redirect output;

							   if (!s.mesh)
								   throw pybind11::value_error("Load the mesh first!");
							   if (!s.mesh->is_volume())
								   throw pybind11::value_error("This function works only on volumetric meshes!");
							   if (s.problem->is_scalar())
								   throw pybind11::value_error("Define a tensor problem!");

							   Eigen::MatrixXd result, stresses, mises;

							   Eigen::MatrixXd v_surf;
							   Eigen::MatrixXi f_surf;
							   const double epsilon = 1e-10;

							   {
								   const polyfem::Mesh3D &mesh3d = *dynamic_cast<polyfem::Mesh3D *>(s.mesh.get());
								   Eigen::MatrixXd points(mesh3d.n_vertices(), 3);
								   Eigen::MatrixXi tets(mesh3d.n_cells(), 4);

								   for (int t = 0; t < mesh3d.n_cells(); ++t)
								   {
									   if (mesh3d.n_cell_vertices(t) != 4)
										   throw pybind11::value_error("Works only with tet meshes!");

									   for (int i = 0; i < 4; ++i)
										   tets(t, i) = mesh3d.cell_vertex(t, i);
								   }

								   for (int p = 0; p < mesh3d.n_vertices(); ++p)
									   points.row(p) << mesh3d.point(p);

								   Eigen::MatrixXi f_surf_tmp, _;
								   igl::boundary_facets(tets, f_surf_tmp);
								   igl::remove_unreferenced(points, f_surf_tmp, v_surf, f_surf, _);
							   }

							   if (apply_displacement)
								   s.interpolate_boundary_tensor_function(v_surf, f_surf, s.sol, s.sol, compute_avg, result, stresses, mises, true);
							   else
								   s.interpolate_boundary_tensor_function(v_surf, f_surf, s.sol, compute_avg, result, stresses, mises, true);

							   return py::make_tuple(v_surf, f_surf, result, stresses, mises);
						   },
						   "returns the traction forces, stresses, and von mises computed on the surface", py::arg("apply_displacement") = bool(false), py::arg("compute_avg") = bool(true))

					   ////////////////////////////////////////////////////////////////////////////////////////////
					   .def(
						   "get_sampled_points_frames", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   assert(!s.solution_frames.empty());

							   std::vector<Eigen::MatrixXd> pts;

							   for (const auto &sol : s.solution_frames)
							   {
								   pts.push_back(sol.points);
							   }

							   return pts;
						   },
						   "returns the points frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def(
						   "get_sampled_connectivity_frames", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   assert(!s.solution_frames.empty());

							   std::vector<Eigen::MatrixXi> tets;

							   for (const auto &sol : s.solution_frames)
								   tets.push_back(sol.connectivity);

							   return tets;
						   },
						   "returns the connectivity frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def(
						   "get_sampled_solution_frames", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   assert(!s.solution_frames.empty());

							   std::vector<Eigen::MatrixXd> fun;

							   for (const auto &sol : s.solution_frames)
							   {
								   fun.push_back(sol.solution);
							   }

							   return fun;
						   },
						   "returns the solution frames for a time dependent problem on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def(
						   "get_sampled_mises_frames", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   assert(!s.solution_frames.empty());

							   std::vector<Eigen::MatrixXd> mises;

							   for (const auto &sol : s.solution_frames)
								   mises.push_back(sol.scalar_value);

							   return mises;
						   },
						   "returns the von mises stresses frames on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def(
						   "get_sampled_mises_avg_frames", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   assert(!s.solution_frames.empty());

							   std::vector<Eigen::MatrixXd> mises;

							   for (const auto &sol : s.solution_frames)
								   mises.push_back(sol.scalar_value_avg);

							   return mises;
						   },
						   "returns the von mises stresses per frame averaged around a vertex on a densly sampled mesh, use 'vismesh_rel_area' to control density")

					   .def(
						   "get_boundary_sidesets", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi faces;
							   Eigen::MatrixXd sidesets;

							   s.get_sidesets(points, faces, sidesets);

							   return py::make_tuple(points, faces, sidesets);
						   },
						   "exports get the boundary sideset, edges in 2d or trangles in 3d")
					   .def(
						   "get_body_ids", [](polyfem::State &s)
						   {
							   //    py::scoped_ostream_redirect output;
							   Eigen::MatrixXd points;
							   Eigen::MatrixXi tets;
							   Eigen::MatrixXi el_id;
							   Eigen::MatrixXd discr;

							   s.build_vis_mesh(points, tets, el_id, discr);

							   Eigen::MatrixXd ids(points.rows(), 1);

							   for (int i = 0; i < points.rows(); ++i)
							   {
								   ids(i) = s.mesh->get_body_id(el_id(i));
							   }

							   return py::make_tuple(points, tets, el_id, ids);
						   },
						   "exports get the body ids")
					   .def(
						   "update_dirichlet_boundary", [](polyfem::State &self, const int id, const py::object &val, const bool isx, const bool isy, const bool isz, const std::string &interp)
						   {
							   using namespace polyfem;
							   // py::scoped_ostream_redirect output;
							   const json json_val = val;
							   if (GenericTensorProblem *prob = dynamic_cast<GenericTensorProblem *>(self.problem.get()))
							   {
								   prob->update_dirichlet_boundary(id, json_val, isx, isy, isz, interp);
							   }
							   else if (GenericScalarProblem *prob = dynamic_cast<GenericScalarProblem *>(self.problem.get()))
							   {
								   prob->update_dirichlet_boundary(id, json_val, interp);
							   }
							   else
							   {
								   throw "Updating BC works only for GenericProblems";
							   }
						   },
						   "updates dirichlet boundary", py::arg("id"), py::arg("val"), py::arg("isx") = bool(true), py::arg("isy") = bool(true), py::arg("isz") = bool(true), py::arg("interp") = std::string(""))
					   .def(
						   "update_neumann_boundary", [](polyfem::State &self, const int id, const py::object &val, const std::string &interp)
						   {
							   using namespace polyfem;
							   // py::scoped_ostream_redirect output;
							   const json json_val = val;

							   if (GenericTensorProblem *prob = dynamic_cast<GenericTensorProblem *>(self.problem.get()))
							   {
								   prob->update_neumann_boundary(id, json_val, interp);
							   }
							   else if (GenericScalarProblem *prob = dynamic_cast<GenericScalarProblem *>(self.problem.get()))
							   {
								   prob->update_neumann_boundary(id, json_val, interp);
							   }
							   else
							   {
								   throw "Updating BC works only for GenericProblems";
							   }
						   },
						   "updates neumann boundary", py::arg("id"), py::arg("val"), py::arg("interp") = std::string(""))
					   .def(
						   "update_pressure_boundary", [](polyfem::State &self, const int id, const py::object &val, const std::string &interp)
						   {
							   using namespace polyfem;
							   // py::scoped_ostream_redirect output;
							   const json json_val = val;

							   if (GenericTensorProblem *prob = dynamic_cast<GenericTensorProblem *>(self.problem.get()))
							   {
								   prob->update_pressure_boundary(id, json_val, interp);
							   }
							   else
							   {
								   throw "Updating BC works only for Tensor GenericProblems";
							   }
						   },
						   "updates pressure boundary", py::arg("id"), py::arg("val"), py::arg("interp") = std::string(""))
					   .def(
						   "render", [rendering_lambda](polyfem::State &self, int width, int height, const Eigen::MatrixXd &camera_positionm, const double camera_fov, const double camera_near, const double camera_far, const bool is_perspective, const Eigen::MatrixXd &lookatm, const Eigen::MatrixXd &upm, const Eigen::MatrixXd &ambient_lightm, const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> &lights, const Eigen::MatrixXd &diffuse_colorm, const Eigen::MatrixXd &specular_colorm, const double specular_exponent)
						   {
							   using namespace Eigen;

							   const int problem_dim = self.problem->is_scalar() ? 1 : self.mesh->dimension();

							   Eigen::MatrixXd tmp = self.boundary_nodes_pos;
							   assert(tmp.rows() * problem_dim == self.sol.size());
							   for (int i = 0; i < self.sol.size(); i += problem_dim)
							   {
								   for (int d = 0; d < problem_dim; ++d)
								   {
									   tmp(i / problem_dim, d) += self.sol(i + d);
								   }
							   }

							   return rendering_lambda(tmp, self.boundary_triangles, width, height, camera_positionm, camera_fov, camera_near, camera_far, is_perspective, lookatm, upm, ambient_lightm, lights, diffuse_colorm, specular_colorm, specular_exponent); 
						   },
						   "renders the scene", py::kw_only(), py::arg("width"), py::arg("height"), py::arg("camera_position") = Eigen::MatrixXd(), py::arg("camera_fov") = double(75), py::arg("camera_near") = double(1), py::arg("camera_far") = double(10), py::arg("is_perspective") = bool(true), py::arg("lookat") = Eigen::MatrixXd(), py::arg("up") = Eigen::MatrixXd(), py::arg("ambient_light") = Eigen::MatrixXd(), py::arg("lights") = std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>>(), py::arg("diffuse_color") = Eigen::MatrixXd(), py::arg("specular_color") = Eigen::MatrixXd(), py::arg("specular_exponent") = double(1))
					   .def(
						   "render_extrinsic", [rendering_lambda](polyfem::State &self, const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>> &vertex_face, int width, int height, const Eigen::MatrixXd &camera_positionm, const double camera_fov, const double camera_near, const double camera_far, const bool is_perspective, const Eigen::MatrixXd &lookatm, const Eigen::MatrixXd &upm, const Eigen::MatrixXd &ambient_lightm, const std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>> &lights, const Eigen::MatrixXd &diffuse_colorm, const Eigen::MatrixXd &specular_colorm, const double specular_exponent) {
							   int v_count = 0;
							   int f_count = 0;
							   for (const auto& vf_pair : vertex_face) {
								   v_count += vf_pair.first.rows();
								   f_count += vf_pair.second.rows();
							   }
							   Eigen::MatrixXd vertices(v_count, 3);
							   Eigen::MatrixXi faces(f_count, 3);
							   v_count = 0;
							   f_count = 0;
							   for (const auto& vf_pair : vertex_face) {
								   vertices.block(v_count, 0, vf_pair.first.rows(), 3) = vf_pair.first;
								   faces.block(f_count, 0, vf_pair.second.rows(), 3) = (vf_pair.second.array() + v_count).matrix();
								   v_count += vf_pair.first.rows();
								   f_count += vf_pair.second.rows();
							   }
							   return rendering_lambda(vertices, faces, width, height, camera_positionm, camera_fov, camera_near, camera_far, is_perspective, lookatm, upm, ambient_lightm, lights, diffuse_colorm, specular_colorm, specular_exponent); 
						   },
						   "renders the extrinsic scene", py::kw_only(), py::arg("vertex_face") = std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXi>>(), py::arg("width"), py::arg("height"), py::arg("camera_position") = Eigen::MatrixXd(), py::arg("camera_fov") = double(75), py::arg("camera_near") = double(1), py::arg("camera_far") = double(10), py::arg("is_perspective") = bool(true), py::arg("lookat") = Eigen::MatrixXd(), py::arg("up") = Eigen::MatrixXd(), py::arg("ambient_light") = Eigen::MatrixXd(), py::arg("lights") = std::vector<std::pair<Eigen::MatrixXd, Eigen::MatrixXd>>(), py::arg("diffuse_color") = Eigen::MatrixXd(), py::arg("specular_color") = Eigen::MatrixXd(), py::arg("specular_exponent") = double(1));

	solver.doc() = "Polyfem solver";

	m.def(
		"polyfem_command", [](const std::string &json_file, const std::string &febio_file, const std::string &mesh, const std::string &problem_name, const std::string &scalar_formulation, const std::string &tensor_formulation, const int n_refs, const bool skip_normalization, const std::string &solver, const int discr_order, const bool p_ref, const bool count_flipped_els, const bool force_linear_geom, const double vis_mesh_res, const bool project_to_psd, const int nl_solver_rhs_steps, const std::string &output, const std::string &vtu, const int log_level, const std::string &log_file, const bool is_quiet, const bool export_material_params)
		{
			json in_args = json({});

			if (!json_file.empty())
			{
				std::ifstream file(json_file);

				if (file.is_open())
					file >> in_args;
				else
					throw "unable to open " + json_file + " file";
				file.close();
			}
			else
			{
				in_args["mesh"] = mesh;
				in_args["force_linear_geometry"] = force_linear_geom;
				in_args["n_refs"] = n_refs;
				if (!problem_name.empty())
					in_args["problem"] = problem_name;
				in_args["normalize_mesh"] = !skip_normalization;
				in_args["project_to_psd"] = project_to_psd;

				if (!scalar_formulation.empty())
					in_args["scalar_formulation"] = scalar_formulation;
				if (!tensor_formulation.empty())
					in_args["tensor_formulation"] = tensor_formulation;
				// in_args["mixed_formulation"] = mixed_formulation;

				in_args["discr_order"] = discr_order;
				// in_args["use_spline"] = use_splines;
				in_args["count_flipped_els"] = count_flipped_els;
				in_args["output"] = output;
				in_args["use_p_ref"] = p_ref;
				// in_args["iso_parametric"] = isoparametric;
				// in_args["serendipity"] = serendipity;

				in_args["nl_solver_rhs_steps"] = nl_solver_rhs_steps;

				if (!vtu.empty())
				{
					in_args["export"]["vis_mesh"] = vtu;
					in_args["export"]["wire_mesh"] = polyfem::StringUtils::replace_ext(vtu, "obj");
				}
				if (!solver.empty())
					in_args["solver_type"] = solver;

				if (vis_mesh_res > 0)
					in_args["vismesh_rel_area"] = vis_mesh_res;

				if (export_material_params)
					in_args["export"]["material_params"] = true;
			}

			polyfem::State state;
			state.init_logger(log_file, log_level, is_quiet);
			state.init(in_args);

			if (!febio_file.empty())
				state.load_febio(febio_file, in_args);
			else
				state.load_mesh();
			state.compute_mesh_stats();

			state.build_basis();

			state.assemble_rhs();
			state.assemble_stiffness_mat();

			state.solve_problem();

			state.compute_errors();

			state.save_json();
			state.export_data();
		},
		"runs the polyfem command, internal usage");

	m.def(
		"solve_febio", [](const std::string &febio_file, const std::string &output_path, const int log_level, const py::kwargs &kwargs)
		{
			if (febio_file.empty())
				throw pybind11::value_error("Specify a febio file!");

			// json in_args = opts.is_none() ? json({}) : json(opts);
			json in_args = json(static_cast<py::dict>(kwargs));

			if (!output_path.empty())
			{
				in_args["export"]["paraview"] = output_path;
				in_args["export"]["wire_mesh"] = polyfem::StringUtils::replace_ext(output_path, "obj");
				in_args["export"]["material_params"] = true;
				in_args["export"]["body_ids"] = true;
				in_args["export"]["contact_forces"] = true;
				in_args["export"]["surface"] = true;
			}

			const int discr_order = in_args.contains("discr_order") ? int(in_args["discr_order"]) : 1;

			if (discr_order == 1 && !in_args.contains("vismesh_rel_area"))
				in_args["vismesh_rel_area"] = 1e10;

			polyfem::State state;
			state.init_logger("", log_level, false);
			state.init(in_args);
			state.load_febio(febio_file, in_args);
			state.compute_mesh_stats();

			state.build_basis();

			state.assemble_rhs();
			state.assemble_stiffness_mat();

			state.solve_problem();

			//state.compute_errors();

			state.save_json();
			state.export_data();
		},
		"runs FEBio", py::arg("febio_file"), py::arg("output_path") = std::string(""), py::arg("log_level") = 2);

	m.def(
		"solve", [](const Eigen::MatrixXd &vertices, const Eigen::MatrixXi &cells, const py::object &sidesets_func, const py::list &sidesets_selection, const py::list &body_selection, const py::list &materials, const std::string &pde, const py::list &diriclet_bc, const py::list &neumann_bc, const py::list &pressure_bc, const py::object &rhs, const bool is_time_dependent, const py::dict &expo, const int log_level, const py::kwargs &kwargs)
		{
			std::string log_file = "";
			const bool is2d = vertices.cols() == 2;

			std::unique_ptr<polyfem::State> res = std::make_unique<polyfem::State>();
			polyfem::State &state = *res;
			state.init_logger(log_file, log_level, false);
			// const int kwargs_size = std::distance(kwargs.begin(), kwargs.end());
			json in_args = json(static_cast<py::dict>(kwargs));

			if (!in_args.contains("normalize_mesh"))
				in_args["normalize_mesh"] = false;

			const bool is_scalar = polyfem::AssemblerUtils::is_scalar(pde);
			in_args["scalar_formulation"] = "";
			in_args["tensor_formulation"] = "";

			if (is_scalar)
				in_args["scalar_formulation"] = pde;
			else
				in_args["tensor_formulation"] = pde;
			in_args["problem"] = is_scalar ? "GenericScalar" : "GenericTensor";

			const int sidesets_selection_size = std::distance(sidesets_selection.begin(), sidesets_selection.end());

			if (!sidesets_func.is_none())
			{
				// const std::function<int(const polyfem::RowVectorNd &)>
				// const std::function<int(const polyfem::RowVectorNd &, bool)>
				// const std::function<int(const std::vector<int> &, bool)>

				//TODO!
				//   if (const auto selection = static_cast<py::dict>(sidesets_func))
				//   {
				// 	  for (auto item : selection)
				// 	  {
				// 		  std::cout << "key=" << std::string(py::str(item.first)) << ", "
				// 					<< "value=" << std::string(py::str(item.second)) << std::endl;
				// 	  }
				//   }
			}
			else if (sidesets_selection_size > 0 && !in_args.contains("boundary_sidesets"))
			{
				json selections = json::array();

				for (const auto &d : sidesets_selection)
				{
					selections.push_back(json(d));
				}

				in_args["boundary_sidesets"] = selections;
			}

			const int body_selection_size = std::distance(body_selection.begin(), body_selection.end());
			if (body_selection_size > 0 && !in_args.contains("body_ids"))
			{
				json selections = json::array();

				for (const auto &d : body_selection)
				{
					selections.push_back(json(d));
				}

				in_args["body_ids"] = selections;
			}

			const int materials_size = std::distance(materials.begin(), materials.end());
			if (materials_size > 1 && !in_args.contains("body_params"))
			{

				json materialss = json::array();

				for (const auto &d : materials)
				{
					materialss.push_back(json(d));
				}

				in_args["body_params"] = materialss;
			}
			else if (materials_size == 1 && !in_args.contains("params"))
			{
				in_args["params"] = json(materials[0]);
			}

			if (!in_args.contains("problem_params"))
				in_args["problem_params"] = {};

			if (!in_args["problem_params"].contains("is_time_dependent"))
				in_args["problem_params"]["is_time_dependent"] = is_time_dependent;

			const int diriclet_bc_size = std::distance(diriclet_bc.begin(), diriclet_bc.end());
			if (diriclet_bc_size > 0 && !in_args["problem_params"].contains("dirichlet_boundary"))
			{
				json bcs = json::array();
				for (const auto &d : diriclet_bc)
				{
					bcs.push_back(json(d));
				}
				in_args["problem_params"]["dirichlet_boundary"] = bcs;
			}
			const int neumann_bc_size = std::distance(neumann_bc.begin(), neumann_bc.end());
			if (neumann_bc_size > 0 && !in_args["problem_params"].contains("neumann_boundary"))
			{
				json bcs = json::array();
				for (const auto &d : neumann_bc)
				{
					bcs.push_back(json(d));
				}
				in_args["problem_params"]["neumann_boundary"] = bcs;
			}
			const int pressure_bc_size = std::distance(pressure_bc.begin(), pressure_bc.end());
			if (pressure_bc_size > 0 && !in_args["problem_params"].contains("pressure_boundary"))
			{
				json bcs = json::array();
				for (const auto &d : pressure_bc)
				{
					bcs.push_back(json(d));
				}
				in_args["problem_params"]["pressure_boundary"] = bcs;
			}

			if (!rhs.is_none() && !in_args["problem_params"].contains("rhs"))
				in_args["problem_params"]["rhs"] = json(rhs);

			json export_json = json(expo);

			if (!in_args.contains("export"))
				in_args["export"] = export_json;

			state.init(in_args);

			if (is2d)
				state.mesh = std::make_unique<polyfem::Mesh2D>();
			else
				state.mesh = std::make_unique<polyfem::Mesh3D>();
			state.mesh->build_from_matrices(vertices, cells);
			state.load_mesh();
			state.compute_mesh_stats();

			state.build_basis();

			state.assemble_rhs();
			state.assemble_stiffness_mat();
			state.solve_problem();

			return res;
		},
		"single solve function", py::kw_only(), py::arg("vertices"), py::arg("cells"), py::arg("sidesets_func") = py::none(), py::arg("sidesets_selection") = py::list(), py::arg("body_selection") = py::list(), py::arg("materials") = py::list(), py::arg("pde") = std::string("LinearElasticity"), py::arg("diriclet_bc") = py::list(), py::arg("neumann_bc") = py::list(), py::arg("pressure_bc") = py::list(), py::arg("rhs") = py::none(), py::arg("is_time_dependent") = bool(false), py::arg("export") = py::dict(), py::arg("log_level") = 2);
}
