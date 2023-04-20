/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MODELS_VIBRONIC_U1_H
#define MODELS_VIBRONIC_U1_H

#ifdef DMRG_VIBRONIC

#include "dmrg/models/measurements/local_at.h"
#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"
#include "dmrg/models/vibrational/VibronicIntegralParser.hpp"

/**
 * @brief Class representing an ab-initio vibronic Hamiltonian.
 * 
 * The Hamiltonian is defined by providing as input the parameters
 * entering the definition of the vibronic Hamiltonian.
 * So far, coupling terms only up to the second-order can be given
 * as input (therefore, only the LVC and the QVC models are supported)
 */

template<class Matrix>
class VibronicModel : public model_impl<Matrix, U1>
{
public:
    // Types definition
    using base = model_impl<Matrix, U1>;
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using value_type = typename Matrix::value_type;
    using pos_t = typename Lattice::pos_t;
    using positions_type = typename std::vector<pos_t>;
    using operators_type = typename std::vector<tag_type>;

    /** @brief Class constructor */
    VibronicModel(const Lattice& lattice_, BaseParameters& parameters_)
      : lat(lattice_), parameters(parameters_), tag_handler(new table_type()), L_(parameters_["L"]),
        n_ele_states_(parameters_["vibronic_nstates"]), n_vib_states_(parameters_["vibronic_nmodes"])
    {
        // Variable definition
        int nMax = parameters_["Nmax"];
        int maxCoupling = chem::getIndexDim(chem::Hamiltonian::Vibronic) / 2; // Two indices are required as the electronic state is also included
        op_t ident_vib_op, ident_ele_op, create_ele_op, destroy_ele_op, count_ele_op;
        op_t position_vib_op, momentum_vib_op;
        // Definition of the physical dimensions.
        // Note that the physical dimension for the vibrations coupled the 0 QC with
        // itself since it cannot generate electronic excitations
        phys.resize(2);
        phys[0].insert(std::make_pair(0, nMax));
        phys[1].insert(std::make_pair(0, 1));
        phys[1].insert(std::make_pair(1, 1));
        // Registering the electronic operators
        ident_ele_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_ele_op.insert_block(Matrix(1, 1, 1), 1, 1);
        create_ele_op.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy_ele_op.insert_block(Matrix(1, 1, 1), 1, 0);
        count_ele_op.insert_block(Matrix(1, 1, 1), 1, 1);
        //
        ident_ele = tag_handler->register_op(ident_ele_op, tag_detail::bosonic);
        create_ele = tag_handler->register_op(create_ele_op, tag_detail::bosonic);
        destroy_ele = tag_handler->register_op(destroy_ele_op, tag_detail::bosonic);
        count_ele = tag_handler->register_op(count_ele_op, tag_detail::bosonic);
        // Registering the vibrational operators
        // Note that, unlike the vibrational case, we add here the 1/sqrt(2)
        // factor, which is not included in the integral file.
        // Moreover, we explicitly write q^2 and p^2 - we don't need to be overly generic,
        // since we don't have so far anharmonic couplings.
        Matrix mpos(nMax, nMax, 0.), mmom(nMax, nMax, 0.), mident(nMax, nMax, 0.);
        // Loads the matrices (code repetition can be in principle avoided with the 
        // "none" vibrational Hamiltonian here)
        mident(0, 0) = 1.;
        for (int n=1; n < nMax; ++n) {
            mident(n, n) = 1.;
            mpos(n-1, n) =  std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mpos(n, n-1) =  std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mmom(n-1, n) =  std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mmom(n, n-1) = -std::sqrt(value_type(n))/std::sqrt(value_type(2.));
        }
        position_vib_op.insert_block(mpos, 0, 0);
        momentum_vib_op.insert_block(mmom, 0, 0);
        ident_vib_op.insert_block(mident, 0, 0);
        ident_vib = tag_handler->register_op(ident_vib_op, tag_detail::bosonic);
        // -- Creates the powers of the position/momentum operator --
        auto powersOfPositions_op = VibrationalHelpers<Matrix, U1>::generatePowersOfPositionOperator(maxCoupling, nMax, ident_vib_op, position_vib_op);
        auto powersOfMomentum_op = VibrationalHelpers<Matrix, U1>::generatePowersOfMomentumOperator(maxCoupling, nMax, ident_vib_op, momentum_vib_op);
        positionPowers.resize(maxCoupling+1);
        momentumPowers.resize(maxCoupling+1);
        positionPowers[0] = ident_vib;
        momentumPowers[0] = ident_vib;
        for (int iOrder = 1; iOrder <= maxCoupling; iOrder++) {
            positionPowers[iOrder] = tag_handler->register_op(powersOfPositions_op[iOrder], tag_detail::bosonic);
            momentumPowers[iOrder] = tag_handler->register_op(powersOfMomentum_op[iOrder], tag_detail::bosonic);
        }
    }

    /** @brief Creates the Hamiltonian terms, i.e. fills the terms_ vector */
    void create_terms() override {
        auto hamiltonianTerms = Vibrational::detail::parseIntegralVibronic<value_type>(parameters, lat);
        for (int iSize = 0; iSize < hamiltonianTerms.first.size(); iSize++) {
            positions_type positions;
            operators_type operators;
            // Electronic contribution
            operators.push_back(create_ele);
            operators.push_back(destroy_ele);
            positions.push_back(hamiltonianTerms.first[iSize][0]);
            positions.push_back(hamiltonianTerms.first[iSize][1]);
            // Skips the first two elements, they are the electronic states
            auto uniqueCoefficients = std::set<int>(hamiltonianTerms.first[iSize].begin()+2,
                                                    hamiltonianTerms.first[iSize].end());
            uniqueCoefficients.erase(0);
            // Manages separately the vertical energy 
            bool isConstantTerm = uniqueCoefficients.size() == 1 && hamiltonianTerms.first[iSize][2] == 0;
            if (!isConstantTerm) {
                for (const auto& iSite: uniqueCoefficients) {
                    // Note that here the index of the 
                    auto vibrationalIndex = lat.get_prop<int>("vibindex", 0, std::abs(iSite)-1);
                    positions.push_back(vibrationalIndex);
                    auto numberOfOccurrences = std::count(hamiltonianTerms.first[iSize].begin()+2,
                                                          hamiltonianTerms.first[iSize].end(), iSite);
                    // Hardcoded so far max 2-body coupling elements
                    assert(numberOfOccurrences > 0 && numberOfOccurrences <= 2);
                    if (iSite < 0)
                        operators.push_back(momentumPowers[numberOfOccurrences]);
                    else if (iSite > 0)
                        operators.push_back(positionPowers[numberOfOccurrences]);
                }
            }
            modelHelper<Matrix, U1>::add_term(positions, operators, hamiltonianTerms.second[iSize], tag_handler, this->terms_, true);
        }
    }

    //TODO ALB to check in which part this goes
    void update(BaseParameters const& p) {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
    }

    // Type physical basis depends on the type of the site (0 --> nuclear, 1 --> electron)
    Index<U1> const& phys_dim(size_t type) const { return phys[type]; }

    /** @brief Getter for the identity matrix */
    tag_type identity_matrix_tag(size_t type) const
    {
        tag_type ret;
        if (type == 0)
            ret = ident_vib;
        else if (type == 1)
            ret = ident_ele;
        else
            throw std::runtime_error("Site type not recognized");
        return ret;
    }

    /** @brief Getter for the filling matrix */
    tag_type filling_matrix_tag(size_t type) const
    {
        tag_type ret;
        if (type == 0)
            ret = ident_vib;
        else if (type == 1)
            ret = ident_ele;
        else
            throw std::runtime_error("Site type not recognized");
        return ret;
    }

    // Total quantum number which must be obtained at the end of the MPS. Should be 1 in all cases.
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const { return 1; }

    /** @brief Getter for the operator associated with a given string */
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "q")
            return positionPowers[1];
        else if (name == "p")
            return momentumPowers[1];
        else if (name == "aplus")
            return create_ele;
        else if (name == "a")
            return destroy_ele;
        else if (name == "id")
            return identity_matrix_tag(type);
        else if (name == "fill")
            return identity_matrix_tag(type);
        else
          throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    /** @brief Getter for the operator table */
    table_ptr operators_table() const { return tag_handler; }

    // Possible measuraments associated to the model
    measurements_type measurements() const
    {
        // Types definition
        using op_vec = std::vector<op_t>;
        using bond_element = std::vector<std::pair<op_vec, bool> >;
        // Variable declaration
        measurements_type meas;
        // Ground state population
        if (parameters.is_set("MEASURE[Population]")) {
            for (std::size_t idx = 0; idx < n_ele_states_; idx++) {
                std::string name = "PopulationState"+std::to_string(idx);
                // Generates vectors for the position operators
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back(idx);
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                identities_local.push_back(this->identity_matrix(0));
                identities_local.push_back(this->identity_matrix(1));
                fillings_local.push_back(this->filling_matrix(0));
                fillings_local.push_back(this->filling_matrix(1));
                // Bonds element (the actual operator involved in the measurement)
                bond_element ops;
                op_vec local_op_vec;
                local_op_vec.push_back(tag_handler->get_op(ident_vib));
                local_op_vec.push_back(tag_handler->get_op(count_ele));
                ops.push_back(std::make_pair(local_op_vec, false));
                meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local,
                                                                      fillings_local, ops));
            }
        }
        return meas;
    }

private:
    const Lattice & lat;
    BaseParameters & parameters;
    std::size_t L_, n_ele_states_, n_vib_states_;
    std::vector< Index<U1> > phys;
    std::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident_ele, count_ele, create_ele, destroy_ele;
    /** Tag for the vibrational identity operator */
    tag_type ident_vib;
    /** Tag for the powers of the position/momentum operators */
    std::vector<tag_type> positionPowers, momentumPowers;
};

#endif // DMRG_VIBRONIC

#endif // MODELS_VIBRONIC_U1_H
