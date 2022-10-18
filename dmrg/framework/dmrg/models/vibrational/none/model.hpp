/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
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

#ifndef MODELS_VIBRATIONAL_NONE_H
#define MODELS_VIBRATIONAL_NONE_H

#ifdef DMRG_VIBRATIONAL

#include <set>
#include <sstream>
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalIntegralParser.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"

/**
 * @brief Class implementing the canonical quantization-based vibrational Hamiltonian.
 *
 * In this model, we use the canonical quantization to map the Born-Oppenheimer
 * vibrational Hamiltonian onto the DMRG lattice.
 * This means that we use the Harmonic Oscillator-based ladder operator to
 * express the momentum/position operators as b^\dagger/b operators.
 */

template<class Matrix>
class WatsonHamiltonian : public model_impl<Matrix, TrivialGroup> {
    // Types definition
    using base = model_impl<Matrix, TrivialGroup>;
    using table_ptr = typename base::table_ptr;
    using table_type = typename base::table_type;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using terms_type = typename std::vector<term_descriptor>;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using positions_type = typename std::vector<typename Lattice::pos_t>;
    using operators_type = typename std::vector<tag_type>;
    using value_type = typename Matrix::value_type;
public:

    /**
     * @brief Class constructor
     * @param lattice object representing the DMRG lattice
     * @param parameters container with the DMRG parameters
     * @param verbose if true, prints information regarding the Hamiltonian terms
     */
    WatsonHamiltonian(const Lattice& lattice, BaseParameters& parameters, bool verbose)
        : lattice_(lattice), parameters_(parameters), tag_handler_(new table_type()), physIndices_(0)
    {
        // Model parameters
        if (parameters_["watson_coordinate_type"] == "cartesian") {
            coordinateType_ = WatsonCoordinateType::CartesianNormalModes;
            maquis::cout << " Coordinate type: Cartesian Normal Modes" << std::endl;
        }
        else if (parameters_["watson_coordinate_type"] == "internal") {
            coordinateType_ = WatsonCoordinateType::InternalNormalModes;
            maquis::cout << " Coordinate type: Internal coordinates-based Normal Modes" << std::endl;
        }
        else {
            throw std::runtime_error("Coordinate type not recognized");
        }
        maxCoupling_ = parameters_["watson_max_coupling"];
        int numModes =  parameters_["L"];
        positionPowers_.resize(numModes);
        momentumPowers_.resize(numModes);
        physIndices_.resize(numModes);
        ident_.resize(numModes);

        std::vector<int> nMaxVec = parameters_["Nmax"].as<std::vector<int> >();
        if (nMaxVec.size()!= numModes && nMaxVec.size()!=1){
            throw std::runtime_error("Nmax needs to be either a single integer or a list with lenght L");
        }

        // Loop over all modes
        for (int mode=0; mode < numModes; mode++){
            int nMax;
            if (nMaxVec.size()==numModes) {
                nMax = nMaxVec[mode];
            } else {
                nMax = nMaxVec[0];
            }
            std::cout << "Nmax for mode " << mode << " is " << nMax << std::endl;
            op_t ident_op, position_op, momentum_op;
            std::vector<op_t> powersOfPositions_op, powersOfMomentum_op;
            TrivialGroup::charge C = TrivialGroup::IdentityCharge;
            int overallDimension = nMax + maxCoupling_;
            
            // Here it's where the "physical" basis is defined
            physIndices_[mode].insert(std::make_pair(C, nMax));
            Matrix mpos(overallDimension, overallDimension, 0.), mmom(overallDimension, overallDimension, 0.);
            Matrix mident(overallDimension, overallDimension, 0.);
            // Loads the matrices
            mident(0, 0) = 1.;
            for (int n = 1; n < overallDimension; n++) {
                mpos(n-1, n) = std::sqrt(value_type(n));
                mpos(n, n-1) = std::sqrt(value_type(n));
                mmom(n-1, n) = std::sqrt(value_type(n));
                mmom(n, n-1) = -std::sqrt(value_type(n));
                mident(n, n) = 1.;
            }
            position_op.insert_block(mpos, C,C);
            momentum_op.insert_block(mmom, C,C);
            ident_op.insert_block(mident, C,C);
            // -- Creates the powers of the position/momentum operator --
            powersOfPositions_op = VibrationalHelpers<Matrix, TrivialGroup>::generatePowersOfPositionOperator(maxCoupling_, nMax, ident_op, position_op);
            powersOfMomentum_op = VibrationalHelpers<Matrix, TrivialGroup>::generatePowersOfMomentumOperator(maxCoupling_, nMax, ident_op, momentum_op);
            // -- Create operator tag table --
            ident_op.resize_block(0, nMax, nMax);
            ident_[mode] = tag_handler_->register_op(ident_op, tag_detail::bosonic);
            positionPowers_[mode].resize(maxCoupling_+1);
            momentumPowers_[mode].resize(maxCoupling_+1);
            positionPowers_[mode][0] = ident_[mode];
            momentumPowers_[mode][0] = ident_[mode];
            for (int iOrder = 1; iOrder <= maxCoupling_; iOrder++) {
                positionPowers_[mode][iOrder] = tag_handler_->register_op(powersOfPositions_op[iOrder], tag_detail::bosonic);
                momentumPowers_[mode][iOrder] = tag_handler_->register_op(powersOfMomentum_op[iOrder], tag_detail::bosonic);
            }
        }
    }

    /** @brief Update the model with the new parameters */
    void update(BaseParameters const &p) {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented or this model.");
    }

    /**
     * @brief Method to load the terms.
     * This method populates the [terms_] member with the Hamiltonian coefficients
     */
    void create_terms() override {
        auto hamiltonianTerms = Vibrational::detail::WatsonIntegralParser<value_type>(parameters_, lattice_, coordinateType_);
        for (const auto& iTerms: hamiltonianTerms) {
            positions_type positions;
            operators_type operators;
            auto termVector = std::vector<int>(iTerms.first.begin(), iTerms.first.end());
            auto newEnd = std::remove(termVector.begin(), termVector.end(), 0);
            auto numberOfNonZeroElements = std::distance(termVector.begin(), newEnd);
            std::stable_sort(termVector.begin(), newEnd, [](const auto& iVal, const auto& jVal) {
                return std::abs(iVal) < std::abs(jVal);
            });
            int outerCounter = 0;
            while (outerCounter < numberOfNonZeroElements) {
                int referenceValue = termVector[outerCounter];
                int mode = abs(referenceValue)-1;
                int innerCounter = 0;
                while (termVector[outerCounter+innerCounter] == referenceValue && innerCounter+outerCounter != numberOfNonZeroElements)
                    innerCounter += 1;
                positions.push_back(abs(referenceValue)-1);
                assert(innerCounter > 0 && innerCounter <= maxCoupling_);
                if (referenceValue < 0)
                    operators.push_back(momentumPowers_[mode][innerCounter]);
                else if (referenceValue > 0)
                    operators.push_back(positionPowers_[mode][innerCounter]);
                outerCounter += innerCounter;
            }
            // Final addition of the terms
            auto coefficient = static_cast<value_type>(iTerms.second);
            modelHelper<Matrix, TrivialGroup>::add_term(positions, operators, coefficient, tag_handler_, this->terms_, true);
        }
    }

    /** @brief Getter for the physical dimension of a given type */
    Index<TrivialGroup> const& phys_dim(size_t type) const { return physIndices_[type]; }

    /** @brief Getter for the identity operator */
    tag_type identity_matrix_tag(size_t type) const { return ident_[type]; }

    /** @brief Getter for the filling operator */
    tag_type filling_matrix_tag(size_t type) const { return identity_matrix_tag(type); }

    /** @brief Gets the quantum number associated with the wfn */
    typename TrivialGroup::charge total_quantum_numbers(BaseParameters& parms) const {
        return typename TrivialGroup::charge();
    }

    /**
     * @brief Gets the operator associated with a given string
     * @param name string describing the operator
     * @param type site type for which the operator is returned
     * @return tag_type tag associated with the requested operator
     */
    tag_type get_operator_tag(const std::string& name, size_t type) const {
        if (name == "id")
            return ident_[type];
        else if (name == "fill")
            return ident_[type];
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    /** @brief Getter for the tag_handler */
    table_ptr operators_table() const {
        return tag_handler_;
    }

    /**
     * @brief Measurement associated with the n-mode Hamiltonian class
     * For now, returns an empty container.
     */
    measurements_type measurements() const {
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        measurements_type meas;
        return meas;
    }

private:

    /** Class member indicating the highest value of the Taylor operator */
    int maxCoupling_;
    /** Ref to the lattice object */
    const Lattice& lattice_;
    /** Parameter container */
    BaseParameters& parameters_;
    /** Physical basis */
    std::vector<Index<TrivialGroup>> physIndices_;
    /** Pointer to the tag_handler */
    std::shared_ptr<TagHandler<Matrix, TrivialGroup> >  tag_handler_;
    /** Tags of the elementary operators */
    std::vector<tag_type> ident_;
    /** Tag for the powers of the position/momentum operators */
    std::vector<std::vector<tag_type>> positionPowers_, momentumPowers_;
    /** Type associated with the vibrational coordinates */
    WatsonCoordinateType coordinateType_;
};

#endif // DMRG_VIBRATIONAL

#endif
