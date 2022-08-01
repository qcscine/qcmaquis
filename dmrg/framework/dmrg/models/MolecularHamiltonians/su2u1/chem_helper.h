/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_CHEM_HELPER_H
#define QC_CHEM_HELPER_H

#include "dmrg/models/MolecularHamiltonians/parse_integrals.h"

namespace chem {
namespace detail {

template <typename M, class S>
class ChemHelperSU2
{
public:
    using value_type = typename M::value_type;
    using term_descriptor = ::term_descriptor<value_type>;
    using pos_t = Lattice::pos_t;

    ChemHelperSU2(BaseParameters & parms, Lattice const & lat, std::shared_ptr<TagHandler<M, S> > tag_handler_)
        : tag_handler(tag_handler_)
    {
        boost::tie(idx_, matrix_elements) = parse_integrals<value_type, S>(parms, lat);
        for (std::size_t m=0; m < matrix_elements.size(); ++m) {
            IndexTuple<S, 4> pos;
            std::copy(idx_.row(m).first, idx_.row(m).second, pos.begin());
            coefficients[pos] = matrix_elements[m];
        }
    
    }
    std::vector<value_type> const & getMatrixElements() const { return matrix_elements; }

    alps::numeric::matrix<Lattice::pos_t> const & getIdx() const { return idx_; }
    int idx(int m, int pos) const {
        return idx_(m,pos);
    }

    void commit_terms(std::vector<term_descriptor> & tagterms) 
    {
        using EightTuple = IndexTuple<S, 8>;
        for (const auto& it: two_terms)
            tagterms.push_back(it.second);
        for (const auto& it: three_terms)
            tagterms.push_back(it.second);
        for (typename std::map<EightTuple, term_descriptor>::const_iterator it = four_terms.begin();
                it != four_terms.end(); ++it)
            tagterms.push_back(it->second);
    }
    
    // Collapse terms with identical operators and different scales into one term
    void add_2term(std::vector<term_descriptor> & tagterms, term_descriptor term)
    {
        IndexTuple<S, 4> id({static_cast<int>(term.position(0)), static_cast<int>(term.position(1)),
                             static_cast<int>(term.operator_tag(0)), static_cast<int>(term.operator_tag(1))});
        if (two_terms.count(id) == 0)
            two_terms[id] = term;
        else
            two_terms[id].coeff += term.coeff;
    }

    void add_3term(std::vector<term_descriptor> & tagterms, term_descriptor term)
    {
        using SixTuple = IndexTuple<S, 6>;
        SixTuple id({static_cast<int>(term.position(0)),
                     static_cast<int>(term.position(1)),
                     static_cast<int>(term.position(2)),
                     static_cast<int>(term.operator_tag(0)),
                     static_cast<int>(term.operator_tag(1)),
                     static_cast<int>(term.operator_tag(2))});
        if (three_terms.count(id) == 0 )
            three_terms[id] = term;
        else
            three_terms[id].coeff += term.coeff;
    }

    void add_4term(std::vector<term_descriptor> & tagterms, term_descriptor term)
    {
        using EightTuple = IndexTuple<S, 8>;
        IndexTuple<S, 4> pos(std::array<int, 4>({term.position(0), term.position(1), term.position(2), term.position(3)}));
        IndexTuple<S, 4> ops(std::array<int, 4>({static_cast<int>(term.operator_tag(0)), static_cast<int>(term.operator_tag(1)),
                                                 static_cast<int>(term.operator_tag(2)), static_cast<int>(term.operator_tag(3))}));
        EightTuple id(pos, ops);
        if (four_terms.count(id) == 0 )
            four_terms[id] = term;
        else
            four_terms[id].coeff += term.coeff;
    }

private:
    std::shared_ptr<TagHandler<M, S> > tag_handler;
    std::vector<value_type> matrix_elements;
    alps::numeric::matrix<Lattice::pos_t> idx_;
    std::map<IndexTuple<S, 4>, value_type> coefficients;
    std::map<IndexTuple<S, 4>, term_descriptor> two_terms;
    std::map<IndexTuple<S, 6>, term_descriptor> three_terms;
    std::map<IndexTuple<S, 8>, term_descriptor> four_terms;
};

}
}

#endif