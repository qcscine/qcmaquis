/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "siteshifter.h"

template<class Matrix, class SymmGroup, class TimeEvolver, class Perturber>
SiteShifter<Matrix,SymmGroup,TimeEvolver,Perturber>::SiteShifter(MPS_type& MPS_input, std::shared_ptr<TimeEvolver> timeEvolver,
                                                                 std::shared_ptr<Perturber> perturber, int init_site, bool backPropagate)
  : MPS_(MPS_input), is_forward_(true), back_propagate_(backPropagate), site_ref_(init_site), zerosite_tensor_(),
    time_evolver_(timeEvolver), perturber_(perturber)
{
    assert (site_ref_ >= 0 && site_ref_ < MPS_.size()) ;
}

template<class Matrix, class SymmGroup, class TimeEvolver, class Perturber>
truncation_results SiteShifter<Matrix, SymmGroup, TimeEvolver, Perturber>::shiftSiteTI(MPO_type const& mpo, boundaries_type& left, boundaries_type& right,
                                                                                      double alpha, double cutoff, std::size_t Mmax, std::size_t Mval,
                                                                                      bool perturb_dm)
{
    // -- Variable declaration --
    truncation_results trunc;
    MPSTen_type res;
    MPSTen_type tmp = MPS_[site_ref_];
    // -- SVD truncation of MPSTensor --
    if (is_forward_) {
        // Simple checks
        assert(site_ref_ < MPS_.size()-1 );
        tmp.make_left_paired();
        perturber_->perturb_mps(tmp, site_ref_, true);
        // Actual calculation
        boost::tie(res, trunc) =
          contraction::Engine<Matrix, Matrix, SymmGroup>::predict_new_state_l2r_sweep(tmp, mpo[site_ref_],
                                                                                      left[site_ref_], right[site_ref_+1],
                                                                                      alpha, cutoff, Mmax, perturb_dm, true);
        res.make_left_paired();
        MPS_[site_ref_].make_left_paired();
        gemm(transpose(conjugate(res.data())), MPS_[site_ref_].data(), zerosite_tensor_);
        left[site_ref_+1] = contr::overlap_mpo_left_step(res, res, left[site_ref_], mpo[site_ref_]);
        MPS_[site_ref_] = res;
        MPS_[site_ref_+1].multiply_from_left(zerosite_tensor_);
    }
    else {
        // Simple check
        assert(site_ref_ > 0) ;
        tmp.make_right_paired() ;
        perturber_->perturb_mps(tmp, site_ref_, false);
        // Actual calculation
        boost::tie(res, trunc)
          = contraction::Engine<Matrix, Matrix, SymmGroup>::predict_new_state_r2l_sweep(tmp, mpo[site_ref_],
                                                                                        left[site_ref_], right[site_ref_+1],
                                                                                        alpha, cutoff, Mmax, perturb_dm, true);
        // MPS perturbation
        res.make_right_paired();
        MPS_[site_ref_].make_right_paired();
        gemm(MPS_[site_ref_].data(), transpose(conjugate(res.data())), zerosite_tensor_);
        right[site_ref_] = contr::overlap_mpo_right_step(res, res, right[site_ref_+1], mpo[site_ref_]);
        MPS_[site_ref_-1].multiply_from_right(zerosite_tensor_);
        MPS_[site_ref_] = res;
    }
    return trunc ;
}

template<class Matrix, class SymmGroup, class TimeEvolver, class Perturber>
truncation_results SiteShifter<Matrix, SymmGroup, TimeEvolver, Perturber>::shiftSiteTD(const MPO_type& mpo, boundaries_type& left, boundaries_type& right,
                                                                                      double alpha, double cutoff, std::size_t Mmax, std::size_t Mval,
                                                                                      bool perturb_dm)
{
    // -- Variable declaration --
    truncation_results trunc;
    MPSTen_type res;
    MPSTen_type tmp = MPS_[site_ref_];
    // -- SVD truncation of MPSTensor --
    if (is_forward_) {
        // Simple checks
        assert(site_ref_ < MPS_.size()-1 );
        tmp.make_left_paired();
        perturber_->perturb_mps(tmp, site_ref_, true);
        // Actual calculation
        boost::tie(res, trunc) =
          contraction::Engine<Matrix, Matrix, SymmGroup>::predict_new_state_l2r_sweep(tmp, mpo[site_ref_],
                                                                                      left[site_ref_], right[site_ref_+1],
                                                                                      alpha, cutoff, Mmax, perturb_dm, true);
        res.make_left_paired();
        MPS_[site_ref_].make_left_paired();
        gemm(transpose(conjugate(res.data())), MPS_[site_ref_].data(), zerosite_tensor_);
        left[site_ref_+1] = contr::overlap_mpo_left_step(res, res, left[site_ref_], mpo[site_ref_]);
        // Final update
        if (back_propagate_)
            perform_propagation(mpo[site_ref_], mpo[site_ref_+1], left[site_ref_+1], right[site_ref_+1]);
        MPS_[site_ref_] = res;
        MPS_[site_ref_+1].multiply_from_left(zerosite_tensor_);
    }
    else {
        // Simple check
        assert(site_ref_ > 0) ;
        tmp.make_right_paired() ;
        perturber_->perturb_mps(tmp, site_ref_, false);
        // Actual calculation
        boost::tie(res, trunc)
          = contraction::Engine<Matrix, Matrix, SymmGroup>::predict_new_state_r2l_sweep(tmp, mpo[site_ref_],
                                                                                        left[site_ref_], right[site_ref_+1],
                                                                                        alpha, cutoff, Mmax, perturb_dm, true);
        // MPS perturbation
        res.make_right_paired();
        MPS_[site_ref_].make_right_paired();
        gemm(MPS_[site_ref_].data(), transpose(conjugate(res.data())), zerosite_tensor_);
        right[site_ref_] = contr::overlap_mpo_right_step(res, res, right[site_ref_+1], mpo[site_ref_]);
        // Final update
        if (back_propagate_)
            perform_propagation(mpo[site_ref_-1], mpo[site_ref_], left[site_ref_], right[site_ref_]);
        MPS_[site_ref_-1].multiply_from_right(zerosite_tensor_);
        MPS_[site_ref_] = res;
    }
    return trunc ;
}

template<class Matrix, class SymmGroup, class TimeEvolver, class Perturber>
void SiteShifter<Matrix, SymmGroup, TimeEvolver, Perturber>::perform_propagation(const MPOTen_type& mpo_ten_left, const MPOTen_type& mpo_ten_right,
                                                                                 const boundary& left, const boundary& right)
{
    if (back_propagate_) {
        std::cout << std::endl ;
        std::cout << "Backward propagating the ZeroSiteTensor object " << std::endl ;
        auto zsp = ZeroSiteProblem<Matrix, SymmGroup>(mpo_ten_left, mpo_ten_right, left, right);
        time_evolver_->evolve(zsp, zerosite_tensor_, false, false);
    }
}
