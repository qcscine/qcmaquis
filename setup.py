import os
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext


def get_virtualenv_path():
    """
    Determines path to install compiled binaries to depending if
    installation is performed with venv or conda/mamba environment active.
    """
    if hasattr(sys, "real_prefix"):
        return sys.prefix

    if hasattr(sys, "base_prefix") and sys.base_prefix != sys.prefix:
        return sys.prefix

    if any(
        keyword in sys.prefix
        for keyword in ["conda", "miniconda", "mamba", "micromamba"]
    ):
        return sys.prefix

    return None


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # Set build type based on DEBUG environment variable
        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # Default build configuration
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            "-DBUILD_SYMMETRIES=TwoU1;TwoU1PG;SU2U1;SU2U1PG",
            "-DBUILD_DMRG_EVOLVE=ON",
            "-DPYTHON_BINDINGS=ON",
            "-DBUILD_TRANSCORRELATED_DMRG=ON",
            "-DENABLE_OMP=ON",
        ]

        # Checks for user-defined "TD" flag via environment variable
        if os.environ.get("TD") == "ON":
            print ("Enabling Compilation of Time-Dependent Vibronic Models")
            cmake_args = [arg for arg in cmake_args if not arg.startswith("-DBUILD_SYMMETRIES=")
            ] # Remove existing -DBUILD_SYMMETRIES 
            cmake_args.append("-DBUILD_SYMMETRIES=U1;NU1;NONE")
            cmake_args.append("-DBUILD_VIBRATIONAL=ON")
            cmake_args.append("-DBUILD_VIBRONIC=ON")
        
        # Vibratinoal build type 
        if os.environ.get("VIB") == "ON":
            print("Enabling Compilation of Vibratinoal Models")
            cmake_args = [arg for arg in cmake_args if not arg.startswith("-DBUILD_SYMMETRIES=")
            ] # Remove existing -DBUILD_SYMMETRIES
            qn = os.environ.get("QN")
            cmake_args.append("-DBUILD_SYMMETRIES=U1;NONE;NU1")
            cmake_args.append("-DBUILD_VIBRATIONAL=ON")
            cmake_args.append("-DCMAKE_CXX_STANDARD=17")
            if qn:
                print("with symmetry number ", qn)
                cmake_args.append(f"-DDMRG_NUMSYMM={qn}")

            

        # Allows to add custom flags using the CMAKE_ARGS environment variable
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        # Build using all cores
        build_args = ["-j"]

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)

        venv = get_virtualenv_path()
        if venv is not None:
            print("Copying compiled binaries to virtual environment.")
            subprocess.check_call(["cp", "qcmaquis", f"{venv}/bin"], cwd=build_temp)


setup(
    name="dmrg",
    version="1.0.0",
    author="ETH Zurich, Laboratory of Physical Chemistry, Reiher Group",
    description="",
    long_description="",
    ext_modules=[CMakeExtension("_dmrg")],
    package_dir={"": "src/python"},
    packages=find_packages(where="src/python"),
    setup_requires=["numpy"],
    install_requires=["numpy"],  # Add any of your other dependencies here
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
