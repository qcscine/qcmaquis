/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef MAQUIS_CPPINTERFACE_H
#define MAQUIS_CPPINTERFACE_H
// Interface layer for Bagel and ChronosQ

#include <complex>
#include <vector>
#include <array>

    typedef std::complex<double> V;
    typedef double dV;
    // Set basic simulation parameters
    // nel -- number of electrons
    // L -- number of sites
    // irrep -- irrep (0-based, for QCMaquis)
    // site_types -- array with irrep # for each site (sized L), NULL if all sites should be of irrep 1
    // conv_thresh -- convergence threshold
    // m -- max. bond dimension
    // nsweeps -- number of sweeps

    void qcmaquis_interface_preinit(int nel, int L, int irrep,
                                    int* site_types,
                                    double conv_thresh, int m,
                                    bool meas_2rdm, bool entropy, bool magnetism,
                                    const std::string& init_type, const std::string& optimization,
                                    int ietl_jcd_maxiter,
                                    double ietl_jcd_tol, double truncation_initial,
                                    double truncation_final, double integral_cutoff,
                                    const std::string& twosite_truncation, const std::string& orb_order);


    // Initializes/updates integrals
    // If integrals were not present before, initializes the interface
    // Otherwise, updates the integrals and re-initializes the model
    //
    // integral_indices: indices are as in FCIDUMP, but flattened
    // integral_values: corresponding values of integrals
    // integral_size: size of integral_value and 1/4 of size of integral_indices
    void qcmaquis_interface_update_integrals(std::vector<std::vector<int>> integral_indices, std::vector<V> integral_values, int integral_size);

    // Set checkpoint names correctly for excited states
    void qcmaquis_interface_set_state(int state);

    // Sets the number of sweeps
    void qcmaquis_interface_set_nsweeps(int nsweeps, int ngrowsweeps, int nmainsweeps);

    void qcmaquis_interface_report(const std::string& reportid);

    // Start a new simulation with stored parameters
    void qcmaquis_interface_reset();

    // Run DMRG optimization
    void qcmaquis_interface_optimize();

    V qcmaquis_interface_get_energy();

    // Get sweep statistics for the last sweep
    void qcmaquis_interface_get_iteration_results(int* nsweeps, std::size_t* m, V* truncated_weight, V* truncated_fraction, V* smallest_ev);

    // Get 1-RDM
    void qcmaquis_interface_get_1rdm(std::vector<std::pair<V,std::array<int,2>>> *rdm1); // BAGEL
    void qcmaquis_interface_get_1rdm(V *rdm1); // CQ
    void qcmaquis_interface_get_1rdm(dV *rdm1){}; // CQ

    // Get 2-RDM
    void qcmaquis_interface_get_2rdm(std::vector<std::pair<V,std::array<int,4>>> *rdm2);// BAGEL
    void qcmaquis_interface_get_2rdm(V *rdm2); // CQ
    void qcmaquis_interface_get_2rdm(dV *rdm2){}; // CQ

#endif
