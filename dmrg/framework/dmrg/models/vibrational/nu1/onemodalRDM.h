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

#ifndef MEASUREMENTS_ONEMODALRDM_H
#define MEASUREMENTS_ONEMODALRDM_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/models/model.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/BaseParameters.h"

#include "utils/io.hpp"
#include <iostream>

namespace measurements {

/** @brief Measurement associated with the one-modal RDM. */
template <class Matrix, int N>
class onemodalRDM : public measurement<Matrix, NU1_template<N>> {
public:
    // Type definition
    using SymmGroup = NU1_template<N>;
    using base = measurement<Matrix, SymmGroup>;
    using op_t = typename model_impl<Matrix, SymmGroup>::op_t;
    using pos_t = typename Lattice::pos_t;
    using tag_type = typename model_impl<Matrix, SymmGroup>::tag_type;
    using positions_type = typename std::vector<pos_t>;
    using operators_type = typename std::vector<tag_type>;
    using impl_type = model_impl<Matrix, SymmGroup>;
    using term_descriptor = typename impl_type::term_descriptor;
    using terms_type = typename std::vector<term_descriptor>;
    using value_type = typename Matrix::value_type;
    using tag_handler = typename std::shared_ptr<TagHandler<Matrix, SymmGroup>>;
    //
    using base::labels_num;
    using base::vector_results;

    /**
     * @brief Class constructor:
     * @param lat_: DMRG lattice
     * @param name_: measurement name (needed by the base class)
     * @param op: vector of the tags associated with the terms that are to be added up to the measurement.
     * @param coeffs: scalar coefficients for the operator.
     * @param tagger: tag_handler associated with the model.
     */
    onemodalRDM(Lattice const&  lat_, std::string name_, std::vector< operators_type > &op,
                std::vector<float_t>& coeffs, const std::shared_ptr<TagHandler<Matrix, SymmGroup>>& tagger,
                const std::vector<tag_type>& ident) : base(name_) , lat(lat_)
    {
        int L = lat_.size();
        for (int iSite = 0; iSite < L; iSite++) {
            int i_type = lat.get_prop<int>("type", iSite);
            terms_type terms_;
            for(int k = 0; k < coeffs.size(); ++k) {
                positions.resize(0);
                operators.resize(0);
                positions.push_back(iSite);
                operators.push_back(op[k][i_type]);
                value_type scaling = static_cast<value_type>(coeffs[k]);
                auto term = modelHelper<Matrix, SymmGroup>::arrange_operators(positions, operators, scaling, tagger);
                term.first.coeff = scaling;
                assert(!term.second);
                terms_.push_back(term.first);
            }
            auto mpoMaker = generate_mpo::TaggedMPOMaker<Matrix, SymmGroup>(lat_, ident, ident, ident, tagger, terms_);
            this->mpoVector.push_back(mpoMaker.create_mpo());
        }
    }

    void evaluate(const MPS<Matrix, SymmGroup>& mps,  boost::optional<reduced_mps<Matrix, SymmGroup> const&> = boost::none) {
        int iSite = 0;
        for (const auto& mpoElement: this->mpoVector) {
            labels_num.push_back({iSite});
            vector_results.push_back(expval(mps, mpoElement));
            iSite++;
        }
        /*
        using contr = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;
        // Creates the left boundary
        mps.canonize(site_);
        auto i = mps[site_].row_dim();
        Boundary<Matrix, SymmGroup> leftBoundary(i, i, 1), leftBoundaryNextSite;
        for (int k = 0; k < leftBoundary[0].n_blocks(); ++k)
            for (int iRow = 0; iRow < num_rows(leftBoundary[0][k]); iRow++)
                for (int iCol = 0; iCol < num_cols(leftBoundary[0][k]); iCol++)
                    leftBoundary[0][k](iRow, iCol) = (iRow == iCol) ? static_cast<value_type>(1.) : static_cast<value_type>(0.);
        // Contraction
        leftBoundaryNextSite = contr::overlap_mpo_left_step(mps[site_], mps[site_], leftBoundary, mpo[site_]);
        assert(leftBoundaryNextSite.aux_dim() == 1);
        this->result = leftBoundaryNextSite.trace();
        */
    }

    /** @brief Cloning method */
    measurement<Matrix, SymmGroup>* do_clone() const { return new onemodalRDM(*this); }

private:
    /* Private members */
    std::vector<MPO<Matrix, SymmGroup>> mpoVector;
    const Lattice& lat;
    positions_type positions;
    operators_type operators;
};

} // namespace measurements

#endif
