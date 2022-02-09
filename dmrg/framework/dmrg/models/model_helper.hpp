/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

#ifndef MAQUIS_DMRG_MODEL_HELPER_HPP
#define MAQUIS_DMRG_MODEL_HELPER_HPP

#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/model.h"

template <class Matrix, class SymmGroup>
class modelHelper {
    // Types definition
    using base = model_impl<Matrix, SymmGroup>;
    using value_type = typename Matrix::value_type;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using terms_type = typename std::vector<term_descriptor>;
    using op_t = typename base::op_t;
    using operators_type = typename std::vector<tag_type>;
    using pos_t = typename Lattice::pos_t;
    using positions_type = typename std::vector<pos_t>;

public:

    /**
     * @brief SQ multiplication routine.
     * 
     * This routine is used to take a list of SQ operators and get the
     * tag corresponding to the product of all SQ operators.
     * Note that operators sitting on the same site are merged together.
     * Note also that zero operators are not added.
     * 
     * @param positions vector with the position where each operator acts.
     * @param operators vector with the operators.
     * @param tag_handler map keeping track of the operator <--> tag association
     * @return
     */
     static std::pair<term_descriptor, bool> arrange_operators(const positions_type& positions, const operators_type& operators,
                                                               value_type& scaling, std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler)
     {
        // Safety check
        assert(positions.size() == operators.size());
        bool FoundZero = false;
        // Types definition
        typedef std::pair<pos_t, tag_type> pos_op_t;
        // Variables definition
        term_descriptor term;
        std::vector<pos_op_t> pos_ops;
        std::transform(positions.begin(), positions.end(), operators.begin(), std::back_inserter(pos_ops),
                       std::make_pair<pos_t const &, tag_type const &>);
        std::stable_sort(pos_ops.begin(), pos_ops.end(), generate_mpo::compare<pos_op_t>);
        // Now that the operators are properly sorted, the formation of the tag can
        // start. Note that get_product_tag returns a new tag if the product does
        // not exists, otherwise returns the existing tag
        for (size_t opnr = 0; opnr < pos_ops.size();) {
            tag_type product = pos_ops[opnr].second;
            size_t range_end = opnr + 1;
            while (range_end < pos_ops.size() && pos_ops[range_end].first == pos_ops[opnr].first) {
                value_type scale = 1.0;
                if (tag_handler->product_is_null(pos_ops[range_end].second, product))
                    FoundZero = true;
                std::tie(product, scale) = tag_handler->get_product_tag(pos_ops[range_end].second, product);
                // At the very end, registers the Hermitian conjugate as well
                registerNewHermitianConjugate(product, tag_handler);
                scaling *= scale;
                range_end++;
            }
            term.push_back(std::make_pair(pos_ops[opnr].first, product));
            opnr = range_end;
        }
        return std::make_pair(term, FoundZero);
    }

    /**
     * @brief Adds a single term to the Hamiltonian object
     * @param positions vector with the position where each operator acts.
     * @param operators vector with the operators.
     * @param coeff Scaling factor for the Hamiltonian
     */
    static void add_term(positions_type const& positions, operators_type const& operators, value_type const& coeff,
                         const std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler, terms_type& terms,
                         bool verbose=false)
    {
        value_type scaling = 1.;
        std::pair<term_descriptor, bool> ret = modelHelper<Matrix, SymmGroup>::arrange_operators(positions, operators, scaling, tag_handler);
        if (!ret.second) {
            auto term = ret.first;
            term.coeff = coeff * scaling;
            if (verbose)
                std::cout << term << std::endl;
            terms.push_back(term);
        }
    }

    /**
     * @brief Register all the operators contained in a given vector
     * @param ops vector with the opeerator to be registered
     * @param kind fermioni/bosoonic
     * @param tag_handler map keeping track of the operator <--> tag association
     * @return std::vector<tag_type> vector of the tags associated with ops
     */
    static std::vector<tag_type> register_all_types(const std::vector<op_t> & ops, tag_detail::operator_kind kind,
                                                    std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler)
    {
        std::vector<tag_type> ret;
        for (std::size_t idx = 0; idx < ops.size(); idx++) {
            std::pair<tag_type, value_type> newtag = tag_handler->checked_register(ops[idx], kind);
            assert( newtag.first < tag_handler->size() );
            assert( std::abs(newtag.second - value_type(1.)) == value_type() );
            ret.push_back(newtag.first);
        }
        return ret;
    }

    /**
     * @brief Register a given pair of Hermitian operators
     * 
     * Note that this call implicitly changes the internal structure of the tag_handler.
     * 
     * @param ops vector with the first set of operators
     * @param hermOps vector with the hermitian conjugates of [ops]
     * @param tag_handler map keeping track of the operator <--> tag association
     * @return std::vector<tag_type> vector of the tags associated with ops
     */
    static void registerHermitianConjugates(const std::vector<tag_type> & ops,
                                            const std::vector<tag_type> & hermOps,
                                            std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler)
    {
        assert(ops.size() == hermOps.size());
        for (int idx = 0; idx < ops.size(); idx++)
            tag_handler->hermitian_pair(ops[idx], hermOps[idx]);
    }

    /**
     * @brief Registers the Hermitian Conjugate of an operator
     * 
     * Note that this method assumes that, if an operator is registered in the tag_handler, also the hermitian conjugate
     * is. The method raises an assertion exception otherwise.
     * 
     * @param tag Tag for which the hermitian conjugate should be registered
     * @param tag_handler std::shared_ptr<TagHandler> pointer carrying the operator table.
     */
    static void registerNewHermitianConjugate(const tag_type& tag, std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler) 
    {
        auto operatorType = (tag_handler->is_fermionic(tag)) ? tag_detail::fermionic : tag_detail::bosonic;
        auto hermitianOperator = tag_handler->get_op(tag);
        hermitianOperator.adjoint_inplace();
        bool alreadyPresent = tag_handler->hasRegistered(hermitianOperator);
        // If already present, it may be either the same or they are already in the table
        if (alreadyPresent) {
            auto newTag = tag_handler->checked_register(hermitianOperator, operatorType);
            assert(newTag.first == tag || tag_handler->herm_conj(tag) == newTag.first);
        }
        else {
            auto newTag = tag_handler->register_op(hermitianOperator, operatorType);
            tag_handler->hermitian_pair(tag, newTag);
        }
    }
};

#endif // MAQUIS_DMRG_MODEL_HELPER_HPP
