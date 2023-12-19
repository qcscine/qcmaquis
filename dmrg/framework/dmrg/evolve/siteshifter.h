/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SITESHIFTER_H
#define SITESHIFTER_H

#include <sstream>
#include <algorithm>
#include <numeric>

#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpotensor.h"
#include "dmrg/mp_tensors/zerositeproblem.h"
#include "dmrg/mp_tensors/siteproblem.h"
#include <dmrg/optimize/ietl_lanczos_solver.h>

template<class Matrix, class SymmGroup, class TimeEvolver, class Perturber>
class SiteShifter
{
public:
    /** Types definition */
    using bm_type = block_matrix<Matrix, SymmGroup>;
    using MPS_type = MPS<Matrix, SymmGroup>;
    using MPO_type = MPO<Matrix, SymmGroup>;
    using MPSTen_type = MPSTensor<Matrix, SymmGroup>;
    using MPOTen_type = MPOTensor<Matrix, SymmGroup>;
    using value_type = typename MPS_type::value_type;
    using scalar_type = typename MPS_type::scalar_type;
    using real_type = double;
    using boundary = Boundary<typename storage::constrained<Matrix>::type, SymmGroup>;
    using boundaries_type = std::vector< boundary >;
    using contr = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;


    /** @brief Same as above, but provides also the time-evolver */
    SiteShifter(MPS_type& MPS_input, std::shared_ptr<TimeEvolver> time_evolver, std::shared_ptr<Perturber> perturber,
                int init_site, bool back_propagate);
    
    /** @brief Getters for the site ref */
    auto getSite() const { return site_ref_; };
    void forward_shift() { site_ref_++; };
    void backward_shift() { site_ref_--; };
    void revert() { is_forward_ = !is_forward_; };

    /**
     * @brief Prints the MPS that is being propagated.
     * @param site_index Reference site
     */
    void print_MPS(int site_index) { std::cout << this->MPS_[site_index] << std::endl ; };

    /** @brief Standard site shifting */
    truncation_results shiftSiteTI(const MPO_type& mpo, boundaries_type& left, boundaries_type& right,
                                   double alpha, double cutoff, std::size_t Mmax, std::size_t Mval, bool perturb_dm);

    /** @brief TD site shifting (inclused the back-propagation of the zero-site tensor) */
    truncation_results shiftSiteTD(const MPO_type& mpo, boundaries_type& left, boundaries_type& right,
                                  double alpha, double cutoff, std::size_t Mmax, std::size_t Mval, bool perturb_dm);

 private:
    // -- PRIVATE METHODS --
    void perform_propagation(const MPOTen_type& mpo_ten_left, const MPOTen_type& mpo_ten_right, const boundary& left,
                             const boundary& right);

    // -- PRIVATE ATTRIBUTES --
    bm_type zerosite_tensor_;
    std::shared_ptr<TimeEvolver> time_evolver_;
    std::shared_ptr<Perturber> perturber_;
    MPS_type& MPS_;
    bool is_forward_, back_propagate_;
    int site_ref_;
};

#include "siteshifter.hpp"

#endif
