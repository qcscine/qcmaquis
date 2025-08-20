import os
import subprocess
import sys
import json
from pathlib import Path

try:
    import tomllib
except ImportError:
    try:
        import tomli as tomllib
    except ImportError:
        tomllib = None

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


def load_build_config():
    """Load build configuration from pyproject.toml"""
    pyproject_path = Path(__file__).parent / "pyproject.toml"
    
    if not pyproject_path.exists():
        raise FileNotFoundError(
            f"pyproject.toml not found at {pyproject_path}. "
            "This file is required for build configuration."
        )
    
    if tomllib is None:
        raise ImportError(
            "tomllib is required to parse pyproject.toml. "
            "For Python < 3.11, install it with: pip install tomli"
        )
    
    with open(pyproject_path, "rb") as f:
        return tomllib.load(f)


def get_build_profile():
    """
    Get the build profile from environment variables.
    Priority: BUILD_PROFILE > legacy env vars > default
    """
    profile_name = "default"
    
    # Check BUILD_PROFILE environment variable
    if os.environ.get("BUILD_PROFILE"):
        profile_name = os.environ.get("BUILD_PROFILE")
        print(f"Using BUILD_PROFILE environment variable: {profile_name}")
        return profile_name
    
    # Legacy environment variable support with deprecation warning
    if os.environ.get("TD") == "ON":
        print("WARNING: TD environment variable is deprecated. Use BUILD_PROFILE=time_dependent instead")
        profile_name = "time_dependent"
    elif os.environ.get("VIB") == "ON":
        print("WARNING: VIB environment variable is deprecated. Use BUILD_PROFILE=vibrational instead")
        profile_name = "vibrational"
    
    return profile_name


class CMakeExtension(Extension):
    def __init__(self, name, sourcedir=""):
        Extension.__init__(self, name, sources=[])
        self.sourcedir = os.path.abspath(sourcedir)


class CMakeBuild(build_ext):
    
    def build_extension(self, ext):
        extdir = os.path.abspath(os.path.dirname(self.get_ext_fullpath(ext.name)))
        if not extdir.endswith(os.path.sep):
            extdir += os.path.sep

        # Load configuration
        config = load_build_config()
        profile_name = get_build_profile()
        
        # Get build profiles
        profiles = config.get("tool", {}).get("dmrg_build", {}).get("profiles", {})
        
        if profile_name not in profiles:
            available = ", ".join(profiles.keys())
            raise ValueError(
                f"Unknown build profile: {profile_name}. "
                f"Available profiles: {available}"
            )
        
        profile = profiles[profile_name]
        
        # Print build configuration
        print("=" * 60)
        print(f"Building DMRG with profile: {profile_name}")
        print("=" * 60)
        print("Configuration:")
        for key, value in profile.items():
            print(f"  {key}: {value}")
        print("  OpenMP: Enabled (always on)")
        print("=" * 60)
        
        # Set build type based on DEBUG environment variable
        debug = int(os.environ.get("DEBUG", 0)) if self.debug is None else self.debug
        cfg = "Debug" if debug else "Release"

        # Build CMake arguments from profile
        cmake_args = [
            f"-DCMAKE_LIBRARY_OUTPUT_DIRECTORY={extdir}",
            f"-DPYTHON_EXECUTABLE={sys.executable}",
            f"-DCMAKE_BUILD_TYPE={cfg}",
            f"-DBUILD_SYMMETRIES={profile['symmetries']}",
            f"-DPYTHON_BINDINGS=ON",
        ]
        
        # Add optional flags based on profile
        if profile.get("evolve", False):
            cmake_args.append("-DBUILD_DMRG_EVOLVE=ON")
        
        if profile.get("transcorrelated", False):
            cmake_args.append("-DBUILD_TRANSCORRELATED_DMRG=ON")
        
        if profile.get("vibrational", False):
            cmake_args.append("-DBUILD_VIBRATIONAL=ON")
        
        if profile.get("vibronic", False):
            cmake_args.append("-DBUILD_VIBRONIC=ON")
        
        # OpenMP is always enabled
        cmake_args.append("-DENABLE_OMP=ON")
        
        if "numsymm" in profile:
            cmake_args.append(f"-DDMRG_NUMSYMM={profile['numsymm']}")
            # Legacy QN environment variable support
            if os.environ.get("QN"):
                print(f"WARNING: QN environment variable is deprecated. Using numsymm={profile['numsymm']} from profile")

        # Validate configuration combinations
        self._validate_configuration(profile, profile_name)
        
        # Allow custom CMAKE_ARGS to override/extend
        if "CMAKE_ARGS" in os.environ:
            custom_args = [item for item in os.environ["CMAKE_ARGS"].split(" ") if item]
            print(f"Adding custom CMAKE_ARGS: {' '.join(custom_args)}")
            cmake_args += custom_args

        # Build using all cores
        build_args = ["-j"]

        build_temp = os.path.join(self.build_temp, ext.name)
        if not os.path.exists(build_temp):
            os.makedirs(build_temp)

        print(f"Running CMake with args: {' '.join(cmake_args)}")
        subprocess.check_call(["cmake", ext.sourcedir] + cmake_args, cwd=build_temp)
        subprocess.check_call(["cmake", "--build", "."] + build_args, cwd=build_temp)

        venv = get_virtualenv_path()
        if venv is not None:
            print("Copying compiled binaries to virtual environment.")
            subprocess.check_call(["cp", "qcmaquis", f"{venv}/bin"], cwd=build_temp)
        
        print("=" * 60)
        print(f"Build completed successfully with profile: {profile_name}")
        print("=" * 60)
    
    def _validate_configuration(self, profile, profile_name):
        """Validate that the configuration is consistent"""
        errors = []
        
        # Check for incompatible symmetry combinations
        symmetries = profile.get("symmetries", "").split(";")
        
        # Vibrational calculations typically shouldn't use electronic symmetries
        if profile.get("vibrational", False):
            electronic_syms = {"TwoU1", "TwoU1PG", "SU2U1", "SU2U1PG"}
            if any(sym in electronic_syms for sym in symmetries):
                errors.append(
                    f"Profile {profile_name} enables vibrational mode but includes "
                    f"electronic symmetries. This may not be intended."
                )
        
        # Time-dependent calculations need evolve enabled
        if profile.get("vibronic", False) and not profile.get("evolve", False):
            errors.append(
                f"Profile {profile_name} enables vibronic mode but evolve is disabled. "
                f"Time-dependent calculations require evolve=true."
            )
        
        if errors:
            print("WARNING: Configuration validation found issues:")
            for error in errors:
                print(f"  - {error}")
            print("Continuing with build despite warnings...")


setup(
    ext_modules=[CMakeExtension("_dmrg")],
    cmdclass={"build_ext": CMakeBuild},
    zip_safe=False,
)
