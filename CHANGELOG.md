# Changelog


## Release 4.0

### Major Features

- Added Transcorrelated DMRG functionality.
- Added Python-bindings for electronic DMRG.
- Added PySCF interface.
- Enhanced SRCAS for electronic DMRG.
- Added functionality to calculate one- and two-modal RDMs in n-mode quantization framework.
- Added $n$-Mode quantized vibrational Hamiltonian.
- Added correlation analysis for vibrational basis functions.

### Minor Changes

- Changed name and location of the executable. It is found here: `build/applications/qcmaquis`.
- Added `simluation_type` keyword for distinguishing types of calculation e.g., `optimize` and `evolve`.
- Removed keyword "ipi_sweeps_per_system", instead the standard keyword "nsweeps" is used.
- Cleaned up and improved output format.

## Release 3.1.4

- Fix missing `#include <map>` in `results_collector.h`

## Release 3.1.3

- Added support for DMRG[IP], the inverse power iteration method applied to MPS wave functions.
- Added support for DMRG[FEAST] for large-scale excited-state DMRG calculations.

## Release 3.1.2

- Added support for vibrational Hamiltonian within the canonical quantization framework.
- Added support for excitonic and vibronic Hamiltonians.

## Release 3.1.1

- Fixed bug raising the "pre-existing term" message when constructing the MPO 
- Enhancement of the overlap calculation in MPS-SI

## Release 3.1.0

- Added support for TD-DMRG with time-independent Hamiltonian

## Release 3.0.6

- Bugfix for compilation with BUILD_PREBO=OFF

## Release 3.0.5

- Added support for PreBO DMRG calculations

## Release 3.0.4
    
- Fixed overlap calculations in MPSSI
- Fixed number of sweeps for DMRGSCF excited states
- Got rid of boost::enable_if and boost::shared_ptr in favour of STL
- Small updates in the Fortran interface and for Mac OS X

## Release 3.0.3

- Fix for Fiedler ordering and excited states
- Faster parallel RDM evaluation
- Fixed some crashes for TDM evaluation
- Added several tests for RDMs
- Ability to read necessary parameters from checkpoint files
