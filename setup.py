import os
import subprocess
import sys

from setuptools import Extension, find_packages, setup
from setuptools.command.build_ext import build_ext

def get_virtualenv_path():
    """Used to work out path to install compiled binaries to."""
    if hasattr(sys, 'real_prefix'):
        return sys.prefix

    if hasattr(sys, 'base_prefix') and sys.base_prefix != sys.prefix:
        return sys.prefix

    if 'conda' in sys.prefix:
        return sys.prefix

    if 'micromamba' in sys.prefix:
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

        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"
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
        if "CMAKE_ARGS" in os.environ:
            cmake_args += [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]

        parallel = 10
        build_args = [f"-j{parallel}"]

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
    packages=find_packages(where="python"),
    setup_requires=["numpy"],
    install_requires=["numpy"],  # Add any of your other dependencies here
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
