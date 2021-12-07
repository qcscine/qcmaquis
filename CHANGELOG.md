# Changelog

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
