/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef MAQUIS_CINTERFACE_H
#define MAQUIS_CINTERFACE_H
// C/Fortran compatible Maquis DMRG interface
// For now, only real matrices and SU2U1PG symmgroup is supported

#include <cstddef>

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
    void qcmaquis_interface_preinit(int nel, int L, int spin, int irrep,
                                 const int* site_types, V conv_thresh, int m, int nsweeps,
                                 const int* sweep_m, int nsweepm, const char* project_name,
                                 bool meas_2rdm);

    // Populate parameters from an existing QCMaquis checkpoint file
    void qcmaquis_interface_preinit_checkpoint(const char* checkpoint_name);

    // Initializes/updates integrals
    // If integrals were not present before, initializes the interface
    // Otherwise, updates the integrals and re-initializes the model
    //
    // integral_indices: indices are as in FCIDUMP, but flattened
    // integral_values: corresponding values of integrals
    // integral_size: size of integral_value and 1/4 of size of integral_indices
    void qcmaquis_interface_update_integrals(const int* integral_indices, const V* integral_values, int integral_size);

    // run preliminary calculation to obtain an MPS as a starting guess from CI-DEAS
    // or the Fiedler ordering
    // if Fiedler ordering is present, return the ordering as a string (starting with 1) in fiedler_order_string (note that its length must be correct!)
    // hf_occupation: Array of HF occupations (as 4,3,2,1) for all states, as flattened row-major array of (L*nstates)
    // For CI-DEAS mandatory, for Fiedler ordering optional
    void qcmaquis_interface_run_starting_guess(int nstates, const char* project_name, bool do_fiedler, bool do_cideas, char* fiedler_order_string, int* hf_occupations);

    // Set checkpoint names correctly for excited states
    void qcmaquis_interface_set_state(int state);

    // Sets the number of sweeps
    void qcmaquis_interface_set_nsweeps(int nsweeps);

    // Set an arbitrary QCMaquis parameter
    void qcmaquis_interface_set_param(const char* key, const char* value);

    // Remove an arbitrary QCMaquis parameter
    void qcmaquis_interface_delete_param(const char* key);

    // Start a new simulation with stored parameters
    void qcmaquis_interface_reset();

    // Run DMRG optimization
    void qcmaquis_interface_optimize();

    // get energy
    double qcmaquis_interface_get_energy();

    // Get sweep statistics for the last sweep
    // nsweeps returns the total number of sweeps made so far
    // returns:
    // m -- maximum bond dimension for the last sweep
    // truncated_weight -- sum of all truncated weights for the last sweep, 0 for single-site optimisation
    // truncated_fraction -- sum of all truncated fractions for the last sweep, 0 for single-site optimisation
    // smallest_ev -- smallest cutoff eigenvalue of the last sweep
    void qcmaquis_interface_get_iteration_results(int* nsweeps, std::size_t* m, V* truncated_weight, V* truncated_fraction, V* smallest_ev);

    // Get 1-RDM
    // size: size of values array, size of indices array = 2*size
    void qcmaquis_interface_get_1rdm(int* indices, V* values, int size);

    // Get spin-1-RDM
    void qcmaquis_interface_get_spdm(int* indices, V* values, int size);

    // Get 2-RDM
    // size: size of values array, size of indices array = 4*size
    void qcmaquis_interface_get_2rdm(int* indices, V* values, int size);

    // Get 3-RDM
    // size: size of values array, size of indices array = 6*size
    void qcmaquis_interface_get_3rdm(int* indices, V* values, int size);

    // Get 4-RDM
    // size: size of values array, size of indices array = 8*size
    void qcmaquis_interface_get_4rdm(int* indices, V* values, int size);

    // Measure 3/4-RDM and save it into HDF5 result file (using project name passed to qcmaquis_interface_preinit)
    // state: state index (starting from 0)
    // Warning: this reinitialises the interface and loads a checkpoint corresponding to the state index
    // so this should be called after wavefunction optimisation

    void qcmaquis_interface_measure_and_save_3rdm(int state);
    void qcmaquis_interface_measure_and_save_4rdm(int state);
    void qcmaquis_interface_measure_and_save_trans3rdm(int state, int bra_state);

    // Measure overlap
    double qcmaquis_interface_get_overlap(const char* filename);

    // Redirect stdout to file filename and restore it back
    void qcmaquis_interface_stdout(const char* filename);
    void qcmaquis_interface_restore_stdout();

    // Obtain number of 3 or 4-RDM elements to be evaluated
    // L: number of orbitals
    // bra_neq_ket: true for transition RDM measurement, otherwise false (only for 3-RDM)
    // slice: integer array of 4 indices corresponding to a 4-RDM slice with four first indices fixed
    // for 3-RDM this has 2 or 3 indices
    // if is not nullptr, number of all 3/4-RDM elements will be evaluated, otherwise only the number of elements in the slice
    int qcmaquis_interface_get_3rdm_elements(int L, bool bra_neq_ket, const int* slice);
    int qcmaquis_interface_get_4rdm_elements(int L, const int* slice);

    enum HIRDM_Template { TEMPLATE_4RDM, TEMPLATE_TRANSITION_3RDM };
    // Prepare QCMaquis input template file for distributed calculation of higher-order RDMs
    // filename: name under which the template should be saved
    // state: state index, beginning with 0
    // tpl: specify if we need to calculate 4-RDM or 3-TDM
    // state_j: second state for 3-TDM measurements, unused for 4-RDM
    void qcmaquis_interface_prepare_hirdm_template(const char* filename, int state, HIRDM_Template tpl, int state_j);
}

#endif