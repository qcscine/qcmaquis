/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODEL_HELPER_HPP
#define MAQUIS_DMRG_MODEL_HELPER_HPP

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
    // +-------------------+
    // | ARRANGE_OPERATORS |
    // +-------------------+
    /**
     * This routine is used to take a list of SQ operators and get the
     * corresponding list of tags. Operators centered on the same center are
     * merged together.
     * @param positions
     * @param operators
     * @param tag_handler
     * @return
     */
     static std::pair<term_descriptor, bool> arrange_operators(const positions_type& positions, const operators_type& operators,
                                                               value_type& scaling, std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler) {
                                                               
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
                if (tag_handler->product_is_null(pos_ops[range_end].second, product)) {
                    FoundZero = true;
                }
                std::tie(product, scale) = tag_handler->get_product_tag(pos_ops[range_end].second, product);
                scaling *= scale;
                range_end++;
            }
            term.push_back(std::make_pair(pos_ops[opnr].first, product));
            opnr = range_end;
        }
        // std::cout << "Overall scaling" << std::endl;
        // std::cout << scaling << std::endl;
        return std::make_pair(term, FoundZero);
    }

    /**
     * @brief Adds a single term to the Hamiltonian object
     * @param positions
     * @param operators
     * @param coeff
     */
    static void add_term(positions_type const& positions, operators_type const& operators, value_type const& coeff,
                         const std::shared_ptr<TagHandler<Matrix, SymmGroup>> tag_handler, terms_type& terms) {
        static int count = 0;
        value_type scaling = 1.;
        std::pair<term_descriptor, bool> ret = modelHelper<Matrix, SymmGroup>::arrange_operators(positions, operators, scaling, tag_handler);
        if (!ret.second) {
            count++;
            auto term = ret.first;
            term.coeff = coeff * scaling;
            terms.push_back(term);
            //if (this->verbose) {
            //    std::cout << term << std::endl;
            //    std::cout << "Operator count = " << count << std::endl;
            //}
        }
    }
};

#endif //MAQUIS_DMRG_MODEL_HELPER_HPP
