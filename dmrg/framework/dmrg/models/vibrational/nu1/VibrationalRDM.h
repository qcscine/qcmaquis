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

/**
 * @brief Enum class for managing the operators in vibrational RDMs.
 *
 * We have three possibilities:
 *
 * 1) AddAll: all possible combinations are added.
 * 2) ExcludeAllSame: exclude terms where all operators sit on the same site.
 * 3) ExcludeSame: exclude terms where there are at least 2 ops sitting on the
 *    same site.
 */
enum class VibrationalRDMModality { AddAll, ExcludeAllSame, ExcludeSame};

/** @brief Measurement associated with a vibrational RDM. */
template <class Matrix, int N>
class VibrationalRDM : public measurement<Matrix, NU1_template<N>> {
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
     * @param excludeSameSite: if true, excludes elements of the RDM sitting on the same site.
     */
    VibrationalRDM(Lattice const&  lat_, std::string name_, const std::vector< std::vector< operators_type >>& op,
                   const std::vector<float_t>& coeffs, const std::shared_ptr<TagHandler<Matrix, SymmGroup>>& tagger,
                   const std::vector<tag_type>& ident, VibrationalRDMModality modality)
        : base(name_) , lat(lat_)
    {
        // Global variables calculation
        int L = lat_.size();
        assert(op.size() == coeffs.size());
        for (const auto& iOp: op)
            assert(iOp.size() == op[0].size());
        int numberOfSQOperators = op[0].size();
        int overallCombinations = std::pow(L, numberOfSQOperators);
        // Generates the operators.
        positions_type positions;
        operators_type operators;
        for (int iTerm = 0; iTerm < overallCombinations; iTerm++) {
            // The representation of iTerm in base L gives the index of each operator.
            auto newIntegerRepresentation = this->to_base(iTerm, L, numberOfSQOperators);
            assert(newIntegerRepresentation.size() <= numberOfSQOperators);
            auto setSize = std::set<int>(newIntegerRepresentation.begin(), newIntegerRepresentation.end()).size();
            bool accept = (modality == VibrationalRDMModality::AddAll) ||
                          (modality == VibrationalRDMModality::ExcludeAllSame && setSize != 1) ||
                          (modality == VibrationalRDMModality::ExcludeSame && setSize == newIntegerRepresentation.size());
            if (accept) {
                // Loop over the terms that compose the operator
                terms_type terms_;
                for(int k = 0; k < coeffs.size(); ++k) {
                    positions.resize(0);
                    operators.resize(0);
                    for (int iOp = 0; iOp < numberOfSQOperators; iOp++) {
                        auto iSite = newIntegerRepresentation[iOp];
                        int i_type = lat.get_prop<int>("type", iSite);
                        positions.push_back(iSite);
                        operators.push_back(op[k][iOp][i_type]);
                    }
                    value_type scaling = static_cast<value_type>(coeffs[k]);
                    auto term = modelHelper<Matrix, SymmGroup>::arrange_operators(positions, operators, scaling, tagger);
                    term.first.coeff = scaling;
                    assert(!term.second);
                    terms_.push_back(term.first);
                }
                auto mpoMaker = generate_mpo::TaggedMPOMaker<Matrix, SymmGroup>(lat_, ident, ident, ident, tagger, terms_);
                this->mpoVector.push_back(mpoMaker.create_mpo());
                // Updates variables of the base class used to retrieve data.
                vector_results.push_back(0);
                labels_num.push_back(positions);
            }
        }
    }

    void evaluate(const MPS<Matrix, SymmGroup>& mps,  boost::optional<reduced_mps<Matrix, SymmGroup> const&> = boost::none) {
        int iSite = 0;
        for (const auto& mpoElement: this->mpoVector) {
            vector_results[iSite] = expval(mps, mpoElement);
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
    measurement<Matrix, SymmGroup>* do_clone() const { return new VibrationalRDM(*this); }

private:

    /** @brief Small helper function to convert a decimal number n into an arbitrary base b */
    std::vector<int> to_base(int n, int base, int overallSize)
    {

        std::vector<int> result(overallSize, 0);
        auto jCont = overallSize-1;
        while (n) {
            result[jCont] = n%base;
            n /= base;
            jCont--;
        }
        return result;
    }

    /* Private members */
    std::vector<MPO<Matrix, SymmGroup>> mpoVector;
    const Lattice& lat;
};

} // namespace measurements

#endif
