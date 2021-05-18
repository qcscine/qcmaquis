//
// Created by robin on 17.05.21.
//

#ifndef MAQUIS_DMRG_MODEL_HELPER_HPP
#define MAQUIS_DMRG_MODEL_HELPER_HPP

#include "dmrg/models/model.h"

template <class Matrix, class SymmGroup>
class modelHelper {
    typedef model_impl<Matrix, SymmGroup> base;
    typedef typename Matrix::value_type value_type;
    // Types definition
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;

    typedef typename base::term_descriptor term_descriptor;
    typedef typename std::vector<term_descriptor> terms_type;
    typedef typename base::op_t op_t;
    typedef typename std::vector<tag_type> operators_type;
    typedef typename base::measurements_type measurements_type;
    typedef typename Lattice::pos_t pos_t;                          // position on the lattice
    typedef typename std::vector<pos_t> positions_type;


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
     static std::pair<term_descriptor, bool> arrange_operators(positions_type const &positions,
                                                               operators_type const &operators,
                                                               value_type &scaling,
                                                               std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler) {
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
                boost::tie(product, scale) = tag_handler->get_product_tag(pos_ops[range_end].second, product);
                scaling *= scale;
                range_end++;
            }
            term.push_back(boost::make_tuple(pos_ops[opnr].first, product));
            opnr = range_end;
        }
        // std::cout << "Overall scaling" << std::endl;
        // std::cout << scaling << std::endl;
        return std::make_pair(term, FoundZero);
    } // arrange_operators




    // +----------+
    // | ADD TERM |
    // +----------+
    /**
     * Adds a single term to the Hamiltonian object
     * @param positions
     * @param operators
     * @param coeff
     */
    static void add_term(positions_type const& positions, operators_type const& operators, value_type const& coeff,
                         const std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler, terms_type& terms) {
        static int count = 0;
        value_type scaling = 1.;
        std::pair<term_descriptor, bool> ret = modelHelper<Matrix,NU1>::arrange_operators(positions, operators, scaling, tag_handler);
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
    } //


};

#endif //MAQUIS_DMRG_MODEL_HELPER_HPP
