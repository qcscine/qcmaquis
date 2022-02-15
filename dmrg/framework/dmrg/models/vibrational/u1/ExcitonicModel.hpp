/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MODELS_EXCITONIC_U1_H
#define MODELS_EXCITONIC_U1_H

#ifdef DMRG_VIBRONIC

#include "dmrg/models/measurements/local_at.h"
#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalModelTraitClass.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"
#include "dmrg/models/vibrational/VibronicIntegralParser.hpp"

/**
 * @brief Class representing the Holstein-Hubbard Hamiltonian.
 * 
 * This Hamiltonian describes an excitonic system composed
 * by N excitons, each one with the ground electronic state
 * and one electronically excited state, and a Coulomb term
 * that couples them via nearest-neighbour couplings.
 * 
 * The Hamiltonian is defined by three components:
 *  1) The Harmonic ground-state PES of each monomer.
 *  2) The LVC component of each excited state.
 *  3) The (vibrationally independent) Coulomb coupling.
 */

template<class Matrix>
class HolsteinHubbardExcitonicHamiltonian : public model_impl<Matrix, U1>
{
public:
    // Types definition
    using base = model_impl<Matrix, U1>;
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using value_type = typename Matrix::value_type;
    using pos_t = typename Lattice::pos_t;

    /** 
     * @brief Model representing an Holstein-Hubbard Hamiltonian
     * In the Holstein-Hubbard Hamiltonian, we have multiple monomers, each one described by
     * Harmonic PESs, where the excited states are modelled with the LVC model.
     * Moreover, off-diagonal coordinate-independent electronic coupling terms are present.
     */
    HolsteinHubbardExcitonicHamiltonian (const Lattice& lat_, BaseParameters & model_) 
        : lat(lat_), model(model_), tag_handler(new table_type()), L_(model["L"]), n_ele_states_(model["vibronic_nstates"]),
          n_vib_states_(model["vibronic_nmodes"]), n_particles_(model["n_excitons"]), phys_indexes(0), J_(0.),
          epsilon_(1.), only_nn_(false)
    {
        // Maximum order of the coupling terms that are supported.
        // For the excitonic Hamiltonian, this will be 
        maxCoupling = VibrationalModelTraitClass<U1>::maximumNumberOfCouplings;
        // Vibronic interaction definition
        J_ = model["J_coupling"].as<value_type>();
        epsilon_ = model["J_excitation"].as<value_type>();
        if (model["J_interaction"] == "nn")
            only_nn_ = true;
        // Variable definition
        std::size_t nMax = model["Nmax"];
        op_t ident_vib_op, ident_ele_op;
        op_t create_ele_op, destroy_ele_op, count_ele_op;
        op_t position_vib_op, momentum_vib_op;
        // Definition of the physical dimensions.
        // First we manage the dimensions for the vibrations, the the ones of the nuclei.
        phys_indexes.resize(2) ;
        phys_indexes[0].insert(std::make_pair(0, nMax));
        phys_indexes[1].insert(std::make_pair(0, 1));
        phys_indexes[1].insert(std::make_pair(1, 1));
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
        Matrix mpos(nMax, nMax, 0.), mmom(nMax, nMax, 0.), mident(nMax, nMax, 0.);
        mident(0,0) = 1.;
        for (int n=1; n < nMax; ++n) {
            mident(n,n) = 1.;
            mpos(n-1,n) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mpos(n,n-1) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mmom(n-1,n) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            mmom(n,n-1) = -std::sqrt(value_type(n))/std::sqrt(value_type(2.));
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
        // == Definition of the Hamiltonian ==
        auto hamiltonianTerms = Vibrational::detail::parseIntegralExcitonic<value_type>(model, lat);
        // == Main loop over the monomers ==
        // We first loop over the number of molecules of the aggregate, and then over the
        // terms entering the vibronic Hamiltonian.
        for (int i_body = 0; i_body < n_particles_; i_body++) {
            std::vector<int> vec_jnk(maxCoupling);
            vec_jnk[0] = i_body;
            // Loop over the Hamiltonian terms
            for (int idx = 0; idx < hamiltonianTerms.first.size(); idx++) {
                std::vector<tag_type> operators;
                std::vector<pos_t> positions;
                for (int op_vib = 0; op_vib < maxCoupling; op_vib++) {
                    if (hamiltonianTerms.first[idx][op_vib] < 0) {
                        operators.push_back(momentumPowers[1]);
                        vec_jnk[1] = -hamiltonianTerms.first[idx][op_vib]-1;
                        positions.push_back(lat.get_prop<int>("vibindex", vec_jnk));
                    }
                    else if (hamiltonianTerms.first[idx][op_vib] > 0) {
                        operators.push_back(positionPowers[1]);
                        vec_jnk[1] = hamiltonianTerms.first[idx][op_vib]-1;
                        positions.push_back(lat.get_prop<int>("vibindex", vec_jnk));
                    }
                }
                // Add the count operator for the specific excited states.
                vec_jnk[1] = 0;
                positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                operators.push_back(count_ele);
                // Builds the term of the Hamiltonian
                modelHelper<Matrix, U1>::add_term(positions, operators, hamiltonianTerms.second[idx], tag_handler, this->terms_, true);
            }
        }
        // Add the J term to the Hamiltonian
        std::vector<int> vec_jnk(2);
        for (int i1_body = 0; i1_body < n_particles_; i1_body++) {
            for (int i2_body = 0; i2_body < n_particles_; i2_body++) {
                if (only_nn_ && (i1_body-i2_body == 1 || i2_body-i1_body == 1) || !only_nn_ && i1_body!=i2_body) {
                    std::vector<tag_type> operators;
                    std::vector<pos_t> positions;
                    vec_jnk[0] = i1_body;
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    vec_jnk[0] = i2_body;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(create_ele);
                    operators.push_back(destroy_ele);
                    modelHelper<Matrix, U1>::add_term(positions, operators, J_, tag_handler, this->terms_);
                }
            }
        }
        // On-site term (for now not used, because we have only one electronic state, but this would be
        // crucial when we have more than a single excited state).
        for (int i1_body = 0; i1_body < n_particles_; i1_body++) {
            vec_jnk[0] = i1_body;
            std::vector<tag_type> operators;
            std::vector<pos_t> positions;
            positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
            operators.push_back(count_ele);
            modelHelper<Matrix, U1>::add_term(positions, operators, epsilon_, tag_handler, this->terms_);
        }
    }

    //TODO ALB to check in which part this goes
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
    }

    /** @brief Getter for the physical basis */
    Index<U1> const& phys_dim(size_t type) const { return phys_indexes[type];
    }
    
    /** @brief Identity matrix getter */
    tag_type identity_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if (type == 0)
            ret = ident_vib;
        else
            ret = ident_ele;
        return ret ;
    }
    
    /** @brief Filling matrix getter */
    tag_type filling_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if (type == 0)
          ret = ident_vib;
        else if (type == 1)
          ret = ident_ele;
        else
          throw std::runtime_error("Site type not recognized") ;
        return ret ;
    }

    /** @brief Identity matrix getter */
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        // ALB Note that here we allow at most 1 particle to be excited.
        return 1;
    }

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
        if (model.is_set("MEASURE[Population]")) {
            for (std::size_t idx = 0; idx < n_particles_; idx++) {
                std::string name = "PopulationState"+std::to_string(idx);
                // Generates vectors for the position operators
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back((n_vib_states_+n_ele_states_)*idx);
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                for (std::size_t idx1 = 0; idx1 < 2; idx1++) {
                    identities_local.push_back(this->identity_matrix(idx1));
                    fillings_local.push_back(this->filling_matrix(idx1));
                }
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
    const Lattice& lat;
    bool only_nn_;
    value_type J_, epsilon_ ;
    BaseParameters& model;
    std::size_t L_, n_ele_states_, n_vib_states_, n_particles_;
    std::vector< Index<U1> > phys_indexes;
    std::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident_vib;
    tag_type ident_ele, count_ele, create_ele, destroy_ele;
    /** Tag for the powers of the position/momentum operators */
    std::vector<tag_type> positionPowers, momentumPowers;
    /** Maximum order of many-body coupling */
    int maxCoupling;
};

#endif // DMRG_VIBRONIC

#endif // MODELS_VIBRONIC_U1_H