from .polyfempy import polyfem_command
import argparse


def polyfem():
    parser = argparse.ArgumentParser()

    parser.add_argument("-j", "--json", type=str,
                        default="", help="Simulation JSON file")

    parser.add_argument("-y", "--yaml", type=str,
                        default="", help="Simulation YAML file")

    parser.add_argument("--max_threads", type=int, default=1,
                        help="Maximum number of threads")

    parser.add_argument("-s", "--strict_validation", action='store_true',
                        help="Enables strict validation of input JSON")

    parser.add_argument("--log_level", type=int, default=1,
                        help="Log level 1 debug 2 info")

    parser.add_argument("-o", "--output_dir", type=str,
                        default="", help="Directory for output files")

    args = parser.parse_args()

    polyfem_command(
        json=args.json,
        yaml=args.yaml,
        log_level=args.log_level,
        strict_validation=args.strict_validation,
        max_threads=args.max_threads,
        output_dir=args.output_dir
        )
