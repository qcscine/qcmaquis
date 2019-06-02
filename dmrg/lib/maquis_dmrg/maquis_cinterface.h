/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*
* This software is part of the ALPS Applications, published under the ALPS
* Application License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
*
* You should have received a copy of the ALPS Application License along with
* the ALPS Applications; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/
#ifndef MAQUIS_CINTERFACE_H
#define MAQUIS_CINTERFACE_H
// C/Fortran compatible Maquis DMRG interface
// For now, only real matrices and SU2U1PG symmgroup is supported


extern "C"
{
    typedef double V;
    // Set basic simulation parameters
    // to be called as initialize_dmrg() from the old interface
    // nel -- number of electrons
    // L -- number of sites
    // spin -- spin
    // irrep -- irrep (0-based, for QCMaquis)
    // site_types -- array with irrep # for each site (sized L), NULL if all sites should be of irrep 1
    // conv_thresh -- convergence threshold
    // m -- max. bond dimension
    // nsweeps -- number of sweeps
    // sweep_m -- if different bond dimensions should be used for different sweeps, provide them here
    // if NULL, it is ignored, otherwise if it's not NULL, m is ignored
    // nsweepm -- size of sweep_m
    // project_name -- prefix for checkpoint files
    // integral_indices: indices are as in FCIDUMP, but flattened
    // integral_values: corresponding values of integrals
    // integral_size: size of integral_value and 1/4 of size of integral_indices
    void qcmaquis_interface_init(int nel, int L, int spin, int irrep,
                                 int* site_types, V conv_thresh, int m, int nsweeps,
                                 int* sweep_m, int nsweepm, char* project_name,
                                 int* integral_indices, V* integral_values, int integral_size, bool meas_2rdm);

    // Set checkpoint names correctly for excited states
    void qcmaquis_interface_set_state(int state);

    // Start a new simulation with stored parameters
    void qcmaquis_interface_reset();

    // Run DMRG optimization
    void qcmaquis_interface_optimize();

    double qcmaquis_interface_get_energy();

    // Get 1-RDM
    // size: size of values array, size of indices array = 2*size
    void qcmaquis_interface_get_1rdm(int* indices, V* values, int size);

    // Get 2-RDM
    // size: size of values array, size of indices array = 4*size
    void qcmaquis_interface_get_2rdm(int* indices, V* values, int size);

}

#endif