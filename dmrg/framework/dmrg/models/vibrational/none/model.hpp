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

#include <sstream>
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalIntegralParser.hpp"

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
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using terms_type = typename std::vector<term_descriptor>;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using pos_t = typename Lattice::pos_t;
    using positions_type = typename std::vector<pos_t>;
    using value_type = typename Matrix::value_type;
    using charge_type = typename TrivialGroup::charge;
public:
    
    /**
     * @brief Class constructor
     * @param lattice_ object representing the DMRG lattice
     * @param parameters_ container with the DMRG parameters
     * @param verbose if true, prints information regarding the Hamiltonian terms
     */
    WatsonHamiltonian(const Lattice& lattice_, BaseParameters& parameters_, bool verbose)
        : lattice(lattice_), parameters(parameters_), tag_handler(new table_type()), physIndices(0)
    {
        // Model parameters
        nMax = parameters_["Nmax"];
        op_t ident_op, create_op, destroy_op, count_op, position_op, momentum_op ;
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        // Here it's where the "physical" basis is defined
        physIndices.insert(std::make_pair(C, nMax));
        Matrix mcount(nMax, nMax), mcreate(nMax, nMax), mdestroy(nMax, nMax), mpos(nMax, nMax),
               mmom(nMax, nMax), mident(nMax, nMax);
        for (int n = 0; n < nMax; n++) {
            for (int j = 0; j < nMax; j++) {
                mcount(n,j)   = 0.;
                mcreate(n,j)  = 0.;
                mdestroy(n,j) = 0.;
                mpos(n,j)     = 0.;
                mmom(n,j)     = 0.;
                mident(n,j)   = 0.;
            }
        }
        mident(0,0) = 1.;
        //create annihilation, creation, position, momentum operators
        for (int n=1; n < nMax; n++) {
            mcount(n,n) = n;
            mident(n,n) = 1.;
            mcreate(n-1,n) = std::sqrt(value_type(n));
            mdestroy(n,n-1) = std::sqrt(value_type(n));
            mpos(n-1,n) = std::sqrt(value_type(n));
            mpos(n,n-1) = std::sqrt(value_type(n));
            mmom(n-1,n) = std::sqrt(value_type(n));
            mmom(n,n-1) = -std::sqrt(value_type(n));
        }
        count_op.insert_block(mcount, C,C);
        create_op.insert_block(mcreate, C,C);
        destroy_op.insert_block(mdestroy, C,C);
        position_op.insert_block(mpos, C,C);
        momentum_op.insert_block(mmom, C,C);
        ident_op.insert_block(mident, C,C);
        // -- Create operator tag table --
        ident = tag_handler->register_op(ident_op, tag_detail::bosonic);
        create = tag_handler->register_op(create_op, tag_detail::bosonic);
        destroy = tag_handler->register_op(destroy_op, tag_detail::bosonic);
        count = tag_handler->register_op(count_op, tag_detail::bosonic);
        position = tag_handler->register_op(position_op, tag_detail::bosonic);
        momentum = tag_handler->register_op(momentum_op, tag_detail::bosonic);
    }

    /** @brief Update the model with the new parameters */
    void update(BaseParameters const &p) {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
    }

    void create_terms() override {
        auto Hamiltonian_term = Vibrational::detail::WatsonIntegralParser<value_type>(parameters, lattice);
        /*
        int hamiltonianSize = Hamiltonian_term.first.size();
        for (int iTerm = 0; iTerm < hamiltonianSize; iTerm++) {
            positions_type positions;
            operators_type operators;
            convertLineToOperators(Hamiltonian_term.first[iTerm], positions, operators);
            modelHelper<Matrix, NU1>::add_term(positions, operators, Hamiltonian_term.second[iTerm], tag_handler, this->terms_);
        }
        */
    }

    /** @brief Getter for the physical dimension of a given type */
    Index<TrivialGroup> const& phys_dim(size_t type) const { return physIndices; }

    /** @brief Getter for the identity operator */
    tag_type identity_matrix_tag(size_t type) const { return ident; }

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
        if (name == "n")
            return count;
        else if (name == "bdag")
            return create;
        else if (name == "b")
            return destroy;
        else if (name == "id")
            return ident;
        else if (name == "fill")
            return ident;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    /** @brief Getter for the tag_handler */
    table_ptr operators_table() const {
        return tag_handler;
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
    const Lattice& lattice;
    int lattice_size, num_modes, nMax;
    BaseParameters& parameters;
    Index<TrivialGroup> physIndices;
    std::shared_ptr<TagHandler<Matrix, TrivialGroup> >  tag_handler;
    tag_type ident, create, destroy, count, position, momentum;
    std::vector<tag_type> positionPowers, momentumPowers;
};

#endif // DMRG_VIBRATIONAL

#endif
