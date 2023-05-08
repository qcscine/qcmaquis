/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MODELS_VIBRATIONAL_NONE_H
#define MODELS_VIBRATIONAL_NONE_H

// #ifdef DMRG_VIBRATIONAL

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
        maquis::cout << std::endl;
        maquis::cout << " == CONSTRUCTING WATSON HAMILTONIAN == " << std::endl;
        maquis::cout << std::endl;
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
        // Set the number of indices which are expected in the FCIDUMP
        // Determines also the maximum many-body coupling degree. Per default read in all integrals that are given
        maxCoupling_ = chem::getIndexDim(chem::Hamiltonian::VibrationalCanonical);
        maxManyBodyCoupling_ = (parameters.is_set("watson_max_coupling")) ? parameters["watson_max_coupling"] : maxCoupling_;
        maxInputManyBodyCoupling_ = (parameters.is_set("watson_max_coupling_input")) ? parameters["watson_max_coupling_input"] : maxCoupling_;
        maquis::cout << " - Maximum many-body coupling order supported: " << maxCoupling_ << std::endl;
        maquis::cout << " - Many-body coupling order expected as input: " << maxInputManyBodyCoupling_ << std::endl;
        maquis::cout << " - Maximum many-body coupling order included in the Hamiltonian " << maxManyBodyCoupling_ << std::endl;
        maquis::cout << std::endl;
        int numModes =  parameters_["L"];
        physIndices_.resize(numModes);
        // Analyzes consistency of nMax parameter
        // nMaxVec = parameters_["Nmax"].as<std::vector<int>>();
        auto nMax = parameters_["Nmax"].as<int>();
        // if (nMaxVec.size() == 1) {
            // auto nMax = nMaxVec[0];
        nMaxVec = std::vector<int>(numModes, nMax);
        // }
        // else {
        //     throw std::runtime_error("Nmax needs to be a single integer");
        // }
        // Loads the physical indices
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        for (int iMode = 0; iMode < numModes; iMode++)
            physIndices_[iMode].insert(std::make_pair(C, nMaxVec[iMode]));
        // Decides how many different dimensions there are
        std::set<int> nMaxUnique(nMaxVec.begin(), nMaxVec.end());
        // Loop over all modes
        for (const auto& nMax: nMaxUnique) {
            op_t ident_op, position_op, momentum_op;
            std::vector<op_t> powersOfPositions_op, powersOfMomentum_op;
            int overallDimension = nMax + maxCoupling_;
            // Here it's where the "physical" basis is defined
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
            powersOfPositions_op = VibrationalHelpers<Matrix, TrivialGroup>::generatePowersOfPositionOperator(maxInputManyBodyCoupling_, nMax, ident_op, position_op);
            powersOfMomentum_op = VibrationalHelpers<Matrix, TrivialGroup>::generatePowersOfMomentumOperator(maxInputManyBodyCoupling_, nMax, ident_op, momentum_op);
            // -- Create operator tag table --
            ident_op.resize_block(0, nMax, nMax);
            ident_[nMax] = tag_handler_->register_op(ident_op, tag_detail::bosonic);
            positionPowers_[nMax].resize(maxInputManyBodyCoupling_+1);
            momentumPowers_[nMax].resize(maxInputManyBodyCoupling_+1);
            positionPowers_[nMax][0] = ident_[nMax];
            momentumPowers_[nMax][0] = ident_[nMax];
            for (int iOrder = 1; iOrder <= maxInputManyBodyCoupling_; iOrder++) {
                positionPowers_[nMax][iOrder] = tag_handler_->register_op(powersOfPositions_op[iOrder], tag_detail::bosonic);
                momentumPowers_[nMax][iOrder] = tag_handler_->register_op(powersOfMomentum_op[iOrder], tag_detail::bosonic);
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
        auto hamiltonianTerms = Vibrational::detail::WatsonIntegralParser<value_type>(parameters_, lattice_, coordinateType_, maxCoupling_,
                                                                                      maxManyBodyCoupling_, maxInputManyBodyCoupling_);
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
                assert(innerCounter > 0);
                if (referenceValue < 0)
                    operators.push_back(momentumPowers_[nMaxVec[mode]][innerCounter]);
                else if (referenceValue > 0)
                    operators.push_back(positionPowers_[nMaxVec[mode]][innerCounter]);
                outerCounter += innerCounter;
            }
            assert(operators.size() == positions.size() && positions.size() <= maxInputManyBodyCoupling_);
            // Final addition of the terms
            auto coefficient = static_cast<value_type>(iTerms.second);
            modelHelper<Matrix, TrivialGroup>::add_term(positions, operators, coefficient, tag_handler_, this->terms_, true);
        }
    }

    /** @brief Getter for the physical dimension of a given type */
    Index<TrivialGroup> const& phys_dim(size_t type) const { return physIndices_[type]; }

    /** @brief Getter for the identity operator */
    tag_type identity_matrix_tag(size_t type) const { return ident_.at(nMaxVec[type]); }

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
            return ident_.at(nMaxVec[type]);
        else if (name == "fill")
            return ident_.at(nMaxVec[type]);
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
    int maxCoupling_, maxManyBodyCoupling_, maxInputManyBodyCoupling_;
    /** Ref to the lattice object */
    const Lattice& lattice_;
    /** Parameter container */
    BaseParameters& parameters_;
    /** Physical basis */
    std::vector<Index<TrivialGroup>> physIndices_;
    /** Pointer to the tag_handler */
    std::shared_ptr<TagHandler<Matrix, TrivialGroup> >  tag_handler_;
    /** Tags of the elementary operators */
    std::unordered_map<int, tag_type> ident_;
    /** Tag for the powers of the position/momentum operators */
    std::unordered_map<int, std::vector<tag_type>> positionPowers_, momentumPowers_;
    /** Type associated with the vibrational coordinates */
    WatsonCoordinateType coordinateType_;
    /** Physical dimension for each site */
    std::vector<int> nMaxVec;
};

// #endif // DMRG_VIBRATIONAL

#endif
