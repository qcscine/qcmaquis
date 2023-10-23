/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef QC_CHEM_HELPER_H
#define QC_CHEM_HELPER_H

#include "dmrg/models/chem/parse_integrals.h"

namespace chem {
namespace detail {

    template <typename M, class S>
    class ChemHelperSU2
    {
    public:
        typedef typename M::value_type value_type;
        typedef ::term_descriptor<value_type> term_descriptor;
        typedef Lattice::pos_t pos_t;
        using InputType = double;

        ChemHelperSU2(BaseParameters & parms, Lattice const & lat, std::shared_ptr<TagHandler<M, S> > tag_handler_)
            : tag_handler(tag_handler_)
        {
            boost::tie(idx_, matrix_elements) = parse_integrals<InputType, S>(parms, lat);

            for (std::size_t m=0; m < matrix_elements.size(); ++m) {
                IndexTuple pos;
                std::copy(idx_.row(m).first, idx_.row(m).second, pos.begin());
                coefficients[pos] = matrix_elements[m];
            }
        }

        const auto& getMatrixElements() const { return matrix_elements; }
        alps::numeric::matrix<Lattice::pos_t> const & getIdx() const { return idx_; }

        int idx(int m, int pos) const {
            return idx_(m,pos);
        }

        void commit_terms(std::vector<term_descriptor> & tagterms) {
            for (typename std::map<IndexTuple, term_descriptor>::const_iterator it = two_terms.begin();
                    it != two_terms.end(); ++it)
                tagterms.push_back(it->second);

            for (typename std::map<SixTuple, term_descriptor>::const_iterator it = three_terms.begin();
                    it != three_terms.end(); ++it)
                tagterms.push_back(it->second);

            for (typename std::map<EightTuple, term_descriptor>::const_iterator it = four_terms.begin();
                    it != four_terms.end(); ++it)
                tagterms.push_back(it->second);
        }

        // Collapse terms with identical operators and different scales into one term

        void add_2term(std::vector<term_descriptor> & tagterms, term_descriptor term)
        {
            IndexTuple id(term.position(0), term.position(1), term.operator_tag(0), term.operator_tag(1));
            if (two_terms.count(id) == 0) {
                two_terms[id] = term;
            }
            else
                two_terms[id].coeff += term.coeff;
        }

        void add_3term(std::vector<term_descriptor> & tagterms, term_descriptor term)
        {
            SixTuple id(term.position(0), term.position(1), term.position(2),
                        term.operator_tag(0), term.operator_tag(1), term.operator_tag(2));
            if (three_terms.count(id) == 0 ) {
                three_terms[id] = term;
            }
            else
                three_terms[id].coeff += term.coeff;
        }

        void add_4term(std::vector<term_descriptor> & tagterms, term_descriptor term)
        {
            IndexTuple pos(term.position(0), term.position(1), term.position(2), term.position(3));
            IndexTuple ops(term.operator_tag(0), term.operator_tag(1), term.operator_tag(2), term.operator_tag(3));
            EightTuple id(pos,ops);
            if (four_terms.count(id) == 0 ) {
                four_terms[id] = term;
            }
            else
                four_terms[id].coeff += term.coeff;
        }

    private:

        std::shared_ptr<TagHandler<M, S> > tag_handler;

        std::vector<InputType> matrix_elements;
        alps::numeric::matrix<Lattice::pos_t> idx_;

        std::map<IndexTuple, value_type> coefficients;

        std::map<IndexTuple, term_descriptor> two_terms;
        std::map<SixTuple, term_descriptor> three_terms;
        std::map<EightTuple, term_descriptor> four_terms;
    };
}
}
#endif
