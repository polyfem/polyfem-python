from .polyfempy import polyfem_command
import argparse


def polyfem():
    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--json", type=str,
                        default="", help="Simulation json file")
    parser.add_argument("-m", "--mesh", type=str, default="", help="Mesh path")
    parser.add_argument("-b", "--febio", type=str,
                        default="", help="FEBio file path")

    parser.add_argument("--n_refs", type=int, default=0,
                        help="Number of refinements")
    parser.add_argument("--not_norm", type=bool, default=True,
                        help="Skips mesh normalization")

    parser.add_argument("--problem", type=str,
                        default="", help="Problem name")
    parser.add_argument("--sform", type=str,
                        default="", help="Scalar formulation")
    parser.add_argument("--tform", type=str,
                        default="", help="Tensor formulation")

    parser.add_argument("--solver", type=str, default="", help="Solver to use")

    parser.add_argument("-q", "-p", type=int, default=1,
                        help="Discretization order")
    parser.add_argument("--p_ref", type=bool,
                        default=False, help="Use p refimenet")
    # parser.add_argument("--spline", use_splines, "Use spline for quad/hex meshes");
    parser.add_argument("--count_flipped_els", type=bool,
                        default=False, help="Count flippsed elements")
    parser.add_argument("--lin_geom", type=bool, default=False,
                        help="Force use linear geometric mapping")
    # parser.add_argument("--isoparametric", isoparametric, "Force use isoparametric basis");
    # parser.add_argument("--serendipity", serendipity, "Use of serendipity elements, only for Q2");
    # parser.add_argument("--stop_after_build_basis", stop_after_build_basis, "Stop after build bases");
    parser.add_argument("--vis_mesh_res", type=float,
                        default=-1.0, help="Vis mesh resolution")
    parser.add_argument("--project_to_psd", type=bool,
                        default=False, help="Project local matrices to psd")
    parser.add_argument("--n_incr_load", type=int, default=-
                        1, help="Number of incremeltal load")

    parser.add_argument("--output", type=str, default="",
                        help="Output json file")
    parser.add_argument("--vtu", type=str, default="", help="Vtu output file")

    parser.add_argument("--quiet", type=bool, default=False,
                        help="Disable cout for logging")
    parser.add_argument("--log_file", type=str,
                        default="", help="Log to a file")
    parser.add_argument("--log_level", type=int, default=1,
                        help="Log level 1 debug 2 info")

    parser.add_argument("--export_material_params", type=bool,
                        default=False, help="Export material parameters")

    args = parser.parse_args()

    polyfem_command(
        args.json,
        args.febio,
        args.mesh,
        args.problem,
        args.sform,
        args.tform,
        args.n_refs,
        args.not_norm,
        args.solver,
        args.q,
        args.p_ref,
        args.count_flipped_els,
        args.lin_geom,
        args.vis_mesh_res,
        args.project_to_psd,
        args.n_incr_load,
        args.output,
        args.vtu,
        args.log_level,
        args.log_file,
        args.quiet,
        args.export_material_params)
