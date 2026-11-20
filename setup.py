"""Build the Python package and its native extension."""

import os
import platform
import re
import subprocess
import sys
import sysconfig
from pathlib import Path

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext
from setuptools.command.build_py import build_py
from setuptools.command.sdist import sdist

ROOT = Path(__file__).resolve().parent
BUILDER_MODULES = {"__init__", "model_builder", "id_relationships", "selection_refs"}


def check_package_sources():
    """Fail before packaging if the generated API or builder is missing."""
    required = [
        "python/polyfempy/generated_api/__init__.py",
        "python/polyfempy/generated_api/generated_api.py",
        "python/polyfempy/generated_api/generated_class.py",
        "python/polyfempy/generated_api/builder_api_manifest.json",
        "python/polyfempy/generated_api/class_tree_manifest.json",
        "python/polyfempy/generated_api/generated_api_manifest.json",
        *[f"python-from-jse/generator/{name}.py" for name in sorted(BUILDER_MODULES)],
    ]
    missing = [name for name in required if not (ROOT / name).is_file()]
    if missing:
        raise RuntimeError(
            "Package sources are missing: " + ", ".join(missing) + ".\n"
            "Initialize the submodules and run: python tools/generate_polyfem_api.py"
        )


class BuildPython(build_py):
    """Ship the generated API and the small runtime used by its builder."""

    def run(self):
        check_package_sources()
        super().run()

    def find_package_modules(self, package, package_dir):
        modules = super().find_package_modules(package, package_dir)
        if package == "polyfempy.generated_api.generator":
            return [module for module in modules if module[1] in BUILDER_MODULES]
        return modules


class SourceDistribution(sdist):
    """Require native sources as well as the generated Python package."""

    def run(self):
        check_package_sources()
        if not (ROOT / "polyfem/CMakeLists.txt").is_file():
            raise RuntimeError(
                "PolyFEM sources are missing. Run: git submodule update --init --recursive"
            )
        super().run()


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        super().__init__(name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def get_output_mapping(self):
        """Keep shared runtimes beside the extension in strict editable installs."""
        mapping = super().get_output_mapping()
        for output, source in list(mapping.items()):
            for pattern in ("*.dll", "*.so", "*.so.*", "*.dylib"):
                for library in Path(source).parent.glob(pattern):
                    if library.is_file():
                        mapping[str(Path(output).parent / library.name)] = str(library)
        return mapping

    def run(self):
        try:
            output = subprocess.check_output(["cmake", "--version"], text=True)
        except OSError as exc:
            raise RuntimeError(
                "CMake 3.25 or later is required to build polyfempy."
            ) from exc
        match = re.search(r"version\s+(\d+)\.(\d+)\.(\d+)", output)
        if match is None or tuple(map(int, match.groups())) < (3, 25, 0):
            raise RuntimeError("CMake 3.25 or later is required to build polyfempy.")
        for ext in self.extensions:
            self.build_extension(ext)

    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        try:
            n_threads = int(os.environ.get("N_THREADS") or min(os.cpu_count() or 1, 2))
        except ValueError as exc:
            raise RuntimeError("N_THREADS must be a positive integer.") from exc
        if n_threads < 1:
            raise RuntimeError("N_THREADS must be a positive integer.")

        python_executable = os.path.abspath(sys.executable)
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={python_executable}",
            f"-DPython_EXECUTABLE={python_executable}",
            f"-DPython3_EXECUTABLE={python_executable}",
            f"-DPython_EXTENSION_MODULE_SUFFIX={sysconfig.get_config_var('EXT_SUFFIX')}",
            f"-DPYTHON_INCLUDE_DIR={sysconfig.get_path('include')}",
            "-DPOLYSOLVE_WITH_SPECTRA=OFF",
            "-DPOLYSOLVE_WITH_AMGCL=OFF",
            "-DPOLYSOLVE_WITH_UMFPACK=OFF",
            "-DPOLYSOLVE_WITH_MPI=OFF",
            "-DPOLYSOLVE_WITH_HYPRE=OFF",
            "-DPOLYFEM_WITH_MISO=OFF",
            "-DPOLYFEM_PORTABLE_BUILD=ON",
            "-DCMAKE_POLICY_VERSION_MINIMUM=3.5",
            "-DUSE_NANOBIND=ON",
            "-DCMAKE_BUILD_TYPE=Release",
        ]
        import nanobind

        cmake_args.append(
            f"-Dnanobind_DIR={(Path(nanobind.__file__).parent / 'cmake').as_posix()}"
        )
        conda_prefix = os.environ.get("CONDA_PREFIX")
        if conda_prefix:
            prefixes = [conda_prefix, os.environ.get("CMAKE_PREFIX_PATH", "")]
            cmake_args.append("-DCMAKE_PREFIX_PATH=" + ";".join(filter(None, prefixes)))

        build_env = os.environ.copy()
        build_workers = n_threads
        if platform.system() == "Darwin":
            cmake_args.extend(
                ["-DPOLYSOLVE_WITH_MKL=OFF", "-DPOLYSOLVE_WITH_ACCELERATE=ON"]
            )
        if platform.system() == "Windows":
            cmake_args.extend(
                [
                    f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY_RELEASE={extdir}",
                    f"-DCMAKE_RUNTIME_OUTPUT_DIRECTORY_RELEASE={extdir}",
                ]
            )
            generator = os.environ.get("CMAKE_GENERATOR", "")
            cache = Path(self.build_temp) / "CMakeCache.txt"
            if not generator and cache.is_file():
                match = re.search(
                    r"^CMAKE_GENERATOR:INTERNAL=(.+)$",
                    cache.read_text(encoding="utf-8"),
                    re.MULTILINE,
                )
                if match:
                    generator = match.group(1).strip()
            if (
                not generator or generator.startswith("Visual Studio")
            ) and sys.maxsize > 2**32:
                cmake_args.extend(["-A", "x64"])
            if not generator or generator.startswith("Visual Studio"):
                # Build one project at a time to keep cl workers within N_THREADS.
                build_workers = 1
                build_env["_CL_"] = (
                    build_env.get("_CL_", "") + f" /MP{n_threads}"
                ).strip()

        build_command = [
            "cmake",
            "--build",
            ".",
            "--config",
            "Release",
            "--target",
            "polyfempy",
            "--parallel",
            str(build_workers),
            "--verbose",
        ]
        Path(self.build_temp).mkdir(parents=True, exist_ok=True)
        for stage, command in (
            ("configuration", ["cmake", ext.sourcedir, *cmake_args]),
            ("build", build_command),
        ):
            try:
                # Stream the original output so a failed build keeps its diagnostics.
                subprocess.check_call(command, cwd=self.build_temp, env=build_env)
            except subprocess.CalledProcessError as exc:
                raise RuntimeError(
                    f"CMake {stage} failed with exit code {exc.returncode}. "
                    f"See the output above. Build directory: {self.build_temp}"
                ) from exc


setup(
    ext_modules=[CMakeExtension("polyfempy.polyfempy")],
    cmdclass={
        "build_ext": CMakeBuild,
        "build_py": BuildPython,
        "sdist": SourceDistribution,
    },
    options={"editable_wheel": {"mode": "strict"}},
    packages=find_packages(where="python", include=["polyfem", "polyfempy", "polyfempy.*"])
    + ["polyfempy.generated_api.generator"],
    package_dir={
        "": "python",
        "polyfempy.generated_api.generator": "python-from-jse/generator",
    },
)
