/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_TWOMODALRDM_H
#define MAQUIS_DMRG_TWOMODALRDM_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/BaseParameters.h"

//namespace measurements {
//
//template <class Matrix, int N>
//class twomodalRDM : public measurement<Matrix, NU1_template<N>> {
//public:
//    using SymmGroup = NU1_template<N>;
//    using base = measurement<Matrix, SymmGroup>;
//    using op_t = typename model_impl<Matrix, SymmGroup>::op_t;
//    using pos_t = typename Lattice::pos_t;
//    using tag_type = typename model_impl<Matrix, SymmGroup>::tag_type;
//    using positions_type = typename std::vector<pos_t>;
//    using operators_type = typename std::vector<tag_type>;
//    using impl_type = model_impl<Matrix, SymmGroup>;
//    using term_descriptor = typename impl_type::term_descriptor;
//    using terms_type = typename std::vector<term_descriptor>;
//    using value_type = typename Matrix::value_type;
//    using tag_handler = typename std::shared_ptr<TagHandler<Matrix, SymmGroup>>;
//
//public:
//    twomodalRDM(Lattice const& lat_, std::string name_, std::vector<term_descriptor> const& terms_,
//                tag_handler const& tagger, operators_type const& ident, int iSite, int jSite) 
//        : base(name_) , lat(lat_), iSite_(std::min(iSite, jSite)), jSite_(std::max(iSite, jSite)) 
//    {
//        MPO<Matrix, SymmGroup> mpo = make_mpo(lat_, tagger, terms_, ident, false);
//        this->mpo = mpo;
//    }
//
//    void evaluate(const MPS<Matrix, SymmGroup>& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> = boost::none) {
//        this-> result = expval(mps, this->mpo);
//        /*
//        using contr = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;
//        // Creates the left boundary
//        mps.canonize(iSite_);
//        auto i = mps[iSite_].row_dim();
//        Boundary<Matrix, SymmGroup> leftBoundary(i, i, 1), leftBoundaryNextSite;
//        for (int k = 0; k < leftBoundary[0].n_blocks(); ++k)
//            for (int iRow = 0; iRow < num_rows(leftBoundary[0][k]); iRow++)
//                for (int iCol = 0; iCol < num_cols(leftBoundary[0][k]); iCol++)
//                    leftBoundary[0][k](iRow, iCol) = (iRow == iCol) ? static_cast<value_type>(1.) : static_cast<value_type>(0.);
//        // Contraction
//        leftBoundaryNextSite = contr::overlap_mpo_left_step(mps[iSite_], mps[iSite_], leftBoundary, mpo[iSite_]);
//        for (int iInter = iSite_+1; iInter < jSite_; iInter++)
//            //leftBoundaryNextSite[0] = contr::overlap_left_step(mps[iInter], mps[iInter], leftBoundaryNextSite[0]);
//            leftBoundaryNextSite = contr::overlap_mpo_left_step(mps[iInter], mps[iInter], leftBoundaryNextSite, mpo[iInter]);
//        auto tmp = contr::overlap_mpo_left_step(mps[jSite_], mps[jSite_], leftBoundaryNextSite, mpo[jSite_]);
//        this->result = tmp.trace();
//        */
//    }
//
//protected:
//    measurement<Matrix, SymmGroup>* do_clone() const {
//        return new twomodalRDM(*this);
//    }
//
//
//private:
//    MPO<Matrix, SymmGroup> mpo;
//    const Lattice& lat;
//    int iSite_, jSite_;
//};
//
//} // namespace measurements

#endif // MAQUIS_DMRG_TWOMODALRDM_H
