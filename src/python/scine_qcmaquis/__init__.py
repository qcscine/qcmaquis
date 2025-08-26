from .maquis_dmrg import QCMaquis

__all__ = ["QCMaquis"]

# Module-level __getattr__ for lazy loading (Python 3.7+)
def __getattr__(name):
    if name == "DMRGSolver":
        try:
            from .pyscf_interface.pyscf_interface import DMRGSolver
            return DMRGSolver
        except ImportError as e:
            raise ImportError(
                f"PySCF is required to use {name}. "
                "Install it with: pip install scine_qcmaquis[pyscf] "
                "or pip install pyscf"
            ) from e
    raise AttributeError(f"module {__name__!r} has no attribute {name!r}")

# For better IDE support and backwards compatibility
def __dir__():
    return __all__ + ["DMRGSolver"]
