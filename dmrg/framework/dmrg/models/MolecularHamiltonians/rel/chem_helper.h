/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
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

#ifndef REL_QC_CHEM_DETAIL_H
#define REL_QC_CHEM_DETAIL_H

#include "dmrg/models/MolecularHamiltonians/parse_integrals.h"

namespace chem { 
namespace detail {

template <typename Matrix, class SymmGroup>
class RelChemHelper
{
public:
    // Coefficients are always complex in relativistic Hamiltonians
    // using value_type = typename Matrix::value_type;
    using  value_type = std::complex<double>;
    using term_descriptor = ::term_descriptor<value_type>;
    using tag_type = typename TagHandler<Matrix, SymmGroup>::tag_type;
    using pos_t = Lattice::pos_t;

    RelChemHelper(BaseParameters& parms, const Lattice& lat_, const std::vector<tag_type>& ident_,
                  const std::vector<tag_type>& fill_, std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_)
        : lat(lat_), ident(ident_), fill(fill_), tag_handler(tag_handler_)
    {
		boost::tie(idx_, matrix_elements) = parse_integrals<value_type, SymmGroup, chem::Hamiltonian::RelativisticElectronic>(parms, lat);
        for (int m = 0; m < matrix_elements.size(); ++m) {
            IndexTuple<SymmGroup, 4> pos;
			std::copy(idx_.row(m).first, idx_.row(m).second, pos.begin());
            coefficients[pos] = matrix_elements[m];
        }
    }

    std::vector<value_type> & getMatrixElements() { return matrix_elements; }

    int idx(int m, int pos) const { return idx_(m,pos); }

    void commit_terms(std::vector<term_descriptor> & tagterms) {
        for (const auto& it : two_terms)
            tagterms.push_back(it.second);
        for (const auto& it : three_terms)
            tagterms.push_back(it.second);
        for (const auto& it : four_terms)
            tagterms.push_back(it.second);
    }

    void add_term(std::vector<term_descriptor> & tagterms,
                  value_type scale, int p1, int p2, std::vector<tag_type> const & op_1, std::vector<tag_type> const & op_2) {

        auto term = TermMaker<Matrix, SymmGroup>::two_term(false, ident, scale, p1, p2, op_1, op_2, tag_handler, lat);
        IndexTuple<SymmGroup, 4> id({p1, p2, static_cast<int>(op_1[lat.get_prop<typename SymmGroup::subcharge>("type", p1)]),
                                             static_cast<int>(op_2[lat.get_prop<typename SymmGroup::subcharge>("type", p2)])});
        if (two_terms.count(id) == 0)
            two_terms[id] = term;
        else
            two_terms[id].coeff += term.coeff;
    }

    void add_term(std::vector<term_descriptor> & tagterms,
                  value_type scale, int s, int p1, int p2,
                  std::vector<tag_type> const & op_i, std::vector<tag_type> const & op_k,
                  std::vector<tag_type> const & op_l, std::vector<tag_type> const & op_j)
    {
        auto term = TermMaker<Matrix, SymmGroup>::three_term(ident, fill, scale, s, p1, p2, op_i, op_k, op_l, op_j, tag_handler, lat);
        IndexTuple<SymmGroup, 6> id({static_cast<int>(term.position(0)), static_cast<int>(term.position(1)),
                                     static_cast<int>(term.position(2)), static_cast<int>(term.operator_tag(0)),
                                     static_cast<int>(term.operator_tag(1)), static_cast<int>(term.operator_tag(2))});
        if (three_terms.count(id) == 0)
            three_terms[id] = term;
        else
            three_terms[id].coeff += term.coeff;
    }

    void add_term(std::vector<term_descriptor> & tagterms, value_type scale,
                  int i, int k, int l, int j,
                  std::vector<tag_type> const & op_i, std::vector<tag_type> const & op_k,
                  std::vector<tag_type> const & op_l, std::vector<tag_type> const & op_j)
    {
		auto term = TermMaker<Matrix, SymmGroup>::four_term(ident, fill, scale, i, k, l, j, op_i, op_k, op_l, op_j, tag_handler, lat);
		if (i<k)
            std::swap(i,k);
		if (j<l)
            std::swap(j,l);
		IndexTuple<SymmGroup, 8> id(IndexTuple<SymmGroup, 4>(std::array<int, 4>{i,k,l,j}),
                                    IndexTuple<SymmGroup, 4>(std::array<int, 4>{static_cast<int>(op_i[lat.get_prop<typename SymmGroup::subcharge>("type",i)]),
                                                                                static_cast<int>(op_k[lat.get_prop<typename SymmGroup::subcharge>("type",k)]),
                                                                                static_cast<int>(op_l[lat.get_prop<typename SymmGroup::subcharge>("type",l)]),
                                                                                static_cast<int>(op_j[lat.get_prop<typename SymmGroup::subcharge>("type",j)])}));
		if (four_terms.count(id) == 0)
			four_terms[id] = term;
		else
			four_terms[id].coeff += term.coeff;
    }

private:
    const std::vector<tag_type>& ident;
    const std::vector<tag_type>& fill;
    std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    Lattice const & lat;
    std::vector<value_type> matrix_elements;
    alps::numeric::matrix<Lattice::pos_t> idx_;
    std::map<IndexTuple<SymmGroup, 4>, value_type> coefficients;
    std::map<IndexTuple<SymmGroup, 8>, term_descriptor> four_terms;
    std::map<IndexTuple<SymmGroup, 6>, term_descriptor> three_terms;
    std::map<IndexTuple<SymmGroup, 4>, term_descriptor> two_terms;
};

} // namespace detail
} // namespace chem

#endif
