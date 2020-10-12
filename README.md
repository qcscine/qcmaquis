# QCMaquis
QCMaquis is an efficient C++11 implementation of the density matrix renormalization group (DMRG) algorithm for quantum chemical Hamiltonians in its matrix product state-matrix product operator (MPS-MPO) formulation. Quantum-chemical operators represented as matrix product operators (MPOs) provide the necessary flexibility to accommodate Abelian and non-Abelian symmetries as well as the implementation of non-relativistic and relativistic quantum chemical Hamiltonians, respectively, in a unified framework. We have implemented the special unitary group of degree 2 (SU(2)) in the MPO representation of the non-relativistic Hamiltonian to ensure spin conservation.

## Current Features:
  - Optimization of spin-adapted SU(2) MPS wave functions with the DMRG algorithm
  - Non-relativistic and scalar-relativistic quantum-chemical Hamiltonians
  - Calculation of excited states
  - A tool set to analyze the MPS wave function and its quantum entanglement
  - One-, two, three- and four-particle reduced density matrices
  - One-, two- and three-particle reduced transition density matrices
  - DMRG-CI and DMRG-SCF interface to the [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas) program package for:
    - DMRG-SCF calculations with and without reaction field (e.g. PCM)
    - State-specific and state-averaged DMRG-SCF calculations
    - Analytic gradients for state-specific DMRG-SCF calculations
    - MPS state interaction (MPS-SI) for the calculation of spin-orbit coupling matrix elements, electronic and magnetic properties
    - state-specific and quasi-degenerate (multi-state) DMRG-NEVPT2 calculations
  - The current release version can be used together with SCINE autoCAS

## Download

We strongly recommend installing and running QCMaquis together with the OpenMOLCAS program package and the OpenMOLCAS-QCMaquis interface, as described below.

## QCMaquis in [OpenMolcas](https://gitlab.com/Molcas/OpenMolcas)

To install QCMaquis with the OpenMOLCAS interface, download OpenMOLCAS from the official OpenMOLCAS GitLab repository and follow its installation guide. To activate QCMaquis (and NEVPT2) support within OpenMolcas, add the following flags to cmake when configuring OpenMOLCAS (see step 3 of the OpenMOLCAS installation guide):

```$ cmake -DDMRG=ON -DNEVPT2=ON <other CMAKE flags> ../```
### Installation Requirements:
- a C++ compiler, optionally a Fortran compiler for the OpenMOLCAS Fortran interface. Only G++ has been tested.
- GNU make
- [CMake](https://cmake.org)
- [Boost](https://boost.org) >= 1.56
- [GSL](https://www.gnu.org/software/gsl/)
- [HDF5](https://www.hdfgroup.org/downloads/hdf5)
- Linear algebra library, such as BLAS/Lapack, [Intel MKL](https://software.intel.com/content/www/us/en/develop/tools/math-kernel-library.html) or [OpenBLAS](https://www.openblas.net/)
- For documentation: PDFLaTeX

## QCMaquis standalone compilation
For a quick start, type the following into the terminal:
```
git clone https://github.com/qcscine/qcmaquis.git qcmaquis
cd qcmaquis/dmrg
mkdir build
cd build
cmake -D<CMAKE_COMPILE_OPTIONS> ../
make
```
where `<CMAKE_COMPILE_OPTIONS>` can be (optionally) one of the following:

### CMake compile options:
- `BUILD_MANUAL`: Compile and install the QCMaquis manual in PDF.
- `BUILD_OPENMOLCAS_INTERFACE`: Build and install the OpenMOLCAS Fortran interface.
- `TESTS`: Compile tests.
- `LAPACK_64_BIT`: Enable if you use a linear algebra library configured with 64-bit integers (ILP64).
- `BLAS_LAPACK_SELECTOR`: Set the vendor of the linear algebra library: `openblas`,`mkl_sequential`, `mkl_parallel`, `veclib` for Accelerate on Mac OS X, `auto` for autodetection and `manual` for setting the linking directories manually are supported. Default is autodetection, which usually does a good job.
- `BUILD_*`: Build (legacy) utility binaries to perform various operations with matrix product state (MPS) checkpoint files.

A full installation guide including package dependencies can be found in Section 3.1 of our manual.

## Documentation

A detailed installation guide and manual for QCMaquis (and its components) can be found [here](https://scine.ethz.ch/static/download/manuals/qcmaquis_manual.pdf). (outdated!)

## Support

For issues and bug reports with QCMaquis in OpenMOLCAS, please open an issue in the [OpenMOLCAS bug tracker](https://gitlab.com/Molcas/OpenMolcas/-/issues) with a QCMaquis tag. For issues and bug reports with the QCMaquis standalone version, please use the [GitHub QCMaquis issue tracker](https://github.com/qcscine/qcmaquis/issues) or send an e-mail to dmrg@phys.chem.ethz.ch.

## References

For reproducibility reasons, please cite, depending on the actual calculations you carried out, one or more of the following papers in publications that present data produced with QCMaquis:

### General reference to the code:
  - S. Keller, M. Dolfi, M. Troyer, M. Reiher, "An efficient matrix product operator representation of the quantum chemical Hamiltonian", J. Chem. Phys., 2015, 143, 244118. [DOI](https://doi.org/10.1063/1.4939000)
### DMRG-NEVPT2:
  - L. Freitag, S. Knecht, C. Angeli, M. Reiher, Multireference Perturbation Theory with Cholesky Decomposition for the Density Matrix Renormalization Group", J. Chem. Theory Comput., 2017, 13, 451. [DOI](https://doi.org/10.1021/acs.jctc.6b00778)
### MPS-SI:
  - S. Knecht, S. Keller, J. Autschbach, M. Reiher, "A Nonorthogonal State-Interaction Approach for Matrix Product State Wave Functions", J. Chem. Theory Comput., 2016, 12, 5881. [DOI](https://doi.org/10.1021/acs.jctc.6b00889)

## QCMaquis and [ALPS](https://alps.comp-phys.org)
QCMaquis builds upon the ALPS MPS project. The ALPS MPS codes implement the DMRG algorithm for variational ground and low-lying excited state search as well as time evolution of arbitrary one- and two-dimensional models in a matrix-product-state representation. They have been developed at ETH Zurich by Michele Dolfi and Bela Bauer in the group of Matthias Troyer with contributions from Sebastian Keller and Alexandr Kosenkov and at the University of Geneva by Timothée Ewart and Adrian Kantian in the group of Thierry Giamarchi. For further information on the ALPS project, please visit https://alps.comp-phys.org and note the original ALPS MPS paper:
- M. Dolfi, B. Bauer, S. Keller, A. Kosenkov, T. Ewart, A. Kantian, T. Giamarchi, M. Troyer, "Matrix product state applications for the ALPS project" ,Comp. Phys. Commun., 2014, 12, 3430. [DOI](https://doi.org/10.1016/j.cpc.2014.08.019)

The current QCMaquis release ships with a modified ALPS library based on ALPS 2.3.0 to reduce compile and runtime dependencies. For more information, please read the README.txt in the `dmrg/alps` subdirectory. The ALPS library is subject to the ALPS Library license, which is found in `dmrg/alps/LICENSE.txt`.
