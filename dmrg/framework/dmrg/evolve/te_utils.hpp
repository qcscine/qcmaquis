/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_TE_UTILS_H
#define APP_TE_UTILS_H

#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"

#include "utils/traits.hpp"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/mpo_manip.h"
#include "dmrg/mp_tensors/compression.h"

#include "dmrg/models/generate_mpo/utils.hpp"

#include <vector>
#include <set>
#include <algorithm>

// Return: Bond terms are allowed to be in the same Hamiltonian object if they do not overlap
//         Site terms are splitted among all Hamiltonian objects using that site
template <typename T>
std::vector<std::vector<term_descriptor<T> > > separate_hamil_terms(std::vector<term_descriptor<T> > const & hamil_terms)
{
    typedef std::map<std::size_t, std::set<std::size_t> > pos_where_t;
    typedef term_descriptor<T> term_t;
    typedef std::map<std::size_t, std::vector<term_t> > pos_terms_t;
    
    std::vector<std::vector<term_t> > ret;
    pos_where_t pos_where;
    pos_terms_t pos_terms;
    
    for (int i=0; i<hamil_terms.size(); ++i)
    {            
        if (hamil_terms[i].size() == 1) {
            pos_terms[hamil_terms[i].position(0)].push_back(hamil_terms[i]);
            continue;
        }
        
        term_t term = hamil_terms[i];
        term.canonical_order(); // TODO: check and fix for fermions!!
        bool used = false;
        for (int n=0; n<ret.size() && !used; ++n)
        {
            bool overlap = false;
            for (int j=0; j<ret[n].size() && !overlap; ++j)
            {
                if ( ret[n][j].site_match(term) ) break;
                overlap = ret[n][j].overlap(term);
            }
            
            if (!overlap) {
            	ret[n].push_back(term);
            	for (int k=0; k<term.size(); ++k)
            		pos_where[term.position(k)].insert(n);
                used = true;
            }
        }
        
        if (!used) {
            ret.push_back( std::vector<term_t>(1, term) );
            
        	for (int k=0; k<term.size(); ++k)
        		pos_where[term.position(k)].insert(ret.size()-1);
        }
    }
    
    // Adding site terms to all Hamiltonians acting on site i
    for (typename pos_terms_t::const_iterator it = pos_terms.begin();
         it != pos_terms.end(); ++it)
    {
        double coeff = 1. / pos_where[it->first].size();
        
        for (typename std::set<std::size_t>::const_iterator it2 = pos_where[it->first].begin();
             it2 != pos_where[it->first].end(); ++it2)
            for (int k=0; k<it->second.size(); ++k)
            {
            	term_t tmp_term = it->second[k];
            	tmp_term.coeff *= coeff;
                ret[*it2].push_back(tmp_term);
            }
    }
    
    // Sorting individual Hamiltonians
    for (int n=0; n<ret.size(); ++n)
        std::sort(ret[n].begin(), ret[n].end());
    
    return ret;
}


template <class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> term2block(typename Matrix::value_type const& scale,
                                           block_matrix<Matrix, SymmGroup> const & op1, Index<SymmGroup> const & phys1_i,
                                           block_matrix<Matrix, SymmGroup> const & op2, Index<SymmGroup> const & phys2_i)
{
    block_matrix<Matrix, SymmGroup> bond_op;
    op_kron(phys1_i, phys2_i, op1, op2, bond_op);
    bond_op *= scale;
    return bond_op;
}

template <class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > hamil_to_blocks(Lattice const& lat, Model<Matrix, SymmGroup> const& model)
{
    std::size_t L = lat.size();
    std::vector<block_matrix<Matrix, SymmGroup> > ret_blocks(L-1);
    
    typename Model<Matrix, SymmGroup>::table_ptr tag_handler = model.operators_table();
    
    typedef term_descriptor<typename Matrix::value_type> term_t;
    std::vector<term_t> const& hamil_terms = model.hamiltonian_terms();
    for (int i=0; i<hamil_terms.size(); ++i)
    {
        term_t term = hamil_terms[i];
        term.canonical_order(); // TODO: check and fix for fermions!!
        std::size_t pos1 = term.position(0);
        int type = lat.get_prop<int>("type", pos1);
        int type_p1, type_m1;
        if (pos1 < L-1) type_p1 = lat.get_prop<int>("type", pos1+1);
        if (pos1 > 0) type_m1 = lat.get_prop<int>("type", pos1-1);
        if (term.size() == 1) {
            if (pos1 != 0 && pos1 != L-1)
                term.coeff /= 2.;
            if (pos1 < L-1)
                ret_blocks[pos1] += term2block(term.coeff,
                                               tag_handler->get_op(term.operator_tag(0)), model.phys_dim(type),
                                               model.identity_matrix(type_p1),            model.phys_dim(type_p1));
            if (pos1 > 0)
                ret_blocks[pos1-1] += term2block(term.coeff,
                                               model.identity_matrix(type_m1),            model.phys_dim(type_m1),
                                               tag_handler->get_op(term.operator_tag(0)), model.phys_dim(type));
        } else if (term.size() == 2) {
            assert( std::abs(term.position(0)-term.position(1)) == 1 );
            ret_blocks[pos1] += term2block(term.coeff,
                                           tag_handler->get_op(term.operator_tag(0)), model.phys_dim(type),
                                           tag_handler->get_op(term.operator_tag(1)), model.phys_dim(type_p1));
        }
    }
    
    return ret_blocks;
}

template <class Matrix, class SymmGroup>
class exp_mpo_maker {
    typedef block_matrix<Matrix, SymmGroup> op_t;
    
    typedef term_descriptor<typename Matrix::value_type> term_t;
    typedef OPTable<Matrix, SymmGroup> op_table_t;
    typedef boost::shared_ptr<op_table_t> op_table_ptr;
    typedef typename op_table_t::tag_type tag_type;
    
    typedef boost::tuple<std::size_t, std::size_t, tag_type, typename Matrix::value_type> pretensor_value;
    typedef std::vector<pretensor_value> pretensor_t;
    typedef std::vector<pretensor_t> prempo_t;
    
public:
    exp_mpo_maker(Index<SymmGroup> phys1_, Index<SymmGroup> phys2_,
                  std::size_t pos1_, std::size_t pos2_, Lattice const& lat_)
    : phys1(phys1_), phys2(phys2_)
    , pos1(pos1_), pos2(pos2_), L(pos2-pos1+1)
    , lattice(lat_)
    , type1(lattice.get_prop<int>("type", pos1))
    , type2(lattice.get_prop<int>("type", pos2))
    , n_boso(0), n_ferm(0)
    { }
    
    void add_term(term_t const& term, Model<Matrix, SymmGroup> const& model)
    {
        if (term.size() > 2)
            throw std::runtime_error("time evolution requires at max two-site bond term.");
        
        op_table_ptr op_table = model.operators_table()->get_operator_table();
        
        /// kron product of operators
        op_t bond_op;
        if (term.size() == 2)
            op_kron(model.phys_dim(type1), model.phys_dim(type2), (*op_table)[term.operator_tag(0)], (*op_table)[term.operator_tag(1)], bond_op);
        else if (term.position(0) == pos1)
            op_kron(model.phys_dim(type1), model.phys_dim(type2), (*op_table)[term.operator_tag(0)], model.identity_matrix(type2), bond_op);
        else if (term.position(0) == pos2)
            op_kron(model.phys_dim(type1), model.phys_dim(type2), model.identity_matrix(type1), (*op_table)[term.operator_tag(0)], bond_op);
        else
            throw std::runtime_error("Operator k not matching any valid position.");
        
        if (term.is_fermionic) {
            fermionic_bond += term.coeff * bond_op;
            n_ferm += 1;
        } else {
            bosonic_bond += term.coeff * bond_op;
            n_boso += 1;
        }
    }
    
    MPO<Matrix, SymmGroup> exp_mpo(typename Matrix::value_type const & alpha,
                                   std::vector<tag_type> const& ident, std::vector<tag_type> const& fill,
                                   op_table_ptr op_table) const
    {
        std::size_t maximum_b = 0;
        prempo_t prempo(L);
        
        if (n_boso > 0) {
            std::vector<op_t> left_ops, right_ops;
            exp_and_split(bosonic_bond, alpha, left_ops, right_ops);
            
            maximum_b = add_to_mpo(prempo, maximum_b, left_ops, right_ops, 1., ident, op_table);
        }
        
        if (n_ferm > 0) {
            /// exp(alpha * op1 \otimes op2)
            {
                std::vector<op_t> left_ops, right_ops;
                exp_and_split(fermionic_bond, alpha, left_ops, right_ops);
                
                maximum_b = add_to_mpo(prempo, maximum_b, left_ops, right_ops, .5, ident, op_table);
                maximum_b = add_to_mpo(prempo, maximum_b, left_ops, right_ops, .5, fill, op_table);
            }
            
            /// exp(-alpha * op1 \otimes op2)
            {
                std::vector<op_t> left_ops, right_ops;
                exp_and_split(fermionic_bond, -alpha, left_ops, right_ops);
                
                maximum_b = add_to_mpo(prempo, maximum_b, left_ops, right_ops, .5, ident, op_table);
                maximum_b = add_to_mpo(prempo, maximum_b, left_ops, right_ops, -.5, fill, op_table);
            }
        }
        
        
        MPO<Matrix, SymmGroup> mpo(L);
        for (size_t p=0; p<L; ++p) {
            size_t nrows, ncols;
            using generate_mpo::rcdim;
            boost::tie(nrows, ncols) = rcdim(prempo[p]);
            MPOTensor<Matrix, SymmGroup> tmp(nrows, ncols, prempo[p], op_table);
            using std::swap;
            swap(mpo[p], tmp);
        }
        
        return mpo;
    }
    
private:
    
    void exp_and_split(op_t const& bond_op, typename Matrix::value_type const & alpha,
                       std::vector<op_t> & left_ops, std::vector<op_t> &  right_ops) const
    {
        op_t bond_exp;
        bond_exp = op_exp_hermitian(phys1*phys2, bond_op, alpha);
        bond_exp = reshape_2site_op(phys1, phys2, bond_exp);
        typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type diag_matrix;
        block_matrix<Matrix, SymmGroup> U, V, left, right;
        block_matrix<diag_matrix, SymmGroup> S, Ssqrt;
        svd(bond_exp, U, V, S);
        Ssqrt = sqrt(S);
        gemm(U, Ssqrt, left);
        gemm(Ssqrt, V, right);
        
#ifdef USE_AMBIENT
        select_proc(ambient::actor_t::common);
        for(std::size_t k = 0; k < S.n_blocks(); ++k){
            ambient::numeric::merge(S[k]);
            ambient::numeric::touch(S[k][0]);
        }
        ambient::sync();
#endif
        for(std::size_t k = 0; k < S.n_blocks(); ++k){
            int keep = std::find_if(S[k].diagonal().first, S[k].diagonal().second, boost::lambda::_1 < 1e-10)-S[k].diagonal().first;
            
            left.resize_block(left.basis().lc(k), left.basis().rc(k),
                              left.basis().ls(k), keep);
            right.resize_block(right.basis().lc(k), right.basis().rc(k),
                               keep, right.basis().rs(k));
        }
        
        left_ops  = reshape_right_to_list(phys1, left);
        right_ops = reshape_left_to_list(phys2, right);
        assert(left_ops.size() == right_ops.size());
    }
    
    std::size_t add_to_mpo(prempo_t & prempo, std::size_t maximum_b,
                           std::vector<op_t> const& left_ops, std::vector<op_t> const& right_ops,
                           typename Matrix::value_type s, std::vector<tag_type> fill,  op_table_ptr op_table) const
    {
        for (std::size_t i=0; i<left_ops.size(); ++i)
        {
            std::size_t b = (maximum_b++);
            
            std::pair<tag_type, typename Matrix::value_type> left_tag, right_tag;
            left_tag = op_table->checked_register(left_ops[i]);
            right_tag = op_table->checked_register(right_ops[i]);
            
            s *= left_tag.second * right_tag.second;
            
            prempo[0].push_back  ( pretensor_value(0, b, left_tag.first,  s ) );
            prempo[L-1].push_back( pretensor_value(b, 0, right_tag.first, 1.) );
            for (std::size_t p=1; p<L-1; ++p)
                prempo[p].push_back( pretensor_value(b, b, fill[lattice.get_prop<int>("prop",pos1+p)], 1.) );
        }
        return maximum_b;
    }
    
    Index<SymmGroup> phys1, phys2;
    std::size_t pos1, pos2, L;
    Lattice lattice;
    int type1, type2;
    std::size_t n_boso, n_ferm;
    op_t bosonic_bond, fermionic_bond;
};

// Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
template <class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_exp_mpo(Lattice const& lat, Model<Matrix, SymmGroup> const& model,
                                    std::vector<term_descriptor<typename Matrix::value_type> > const& hamil_terms,
                                    typename Matrix::value_type const & alpha = 1)
{
    typedef term_descriptor<typename Matrix::value_type> term_t;
    typedef OPTable<Matrix, SymmGroup> op_table_t;
    typedef boost::shared_ptr<op_table_t> op_table_ptr;
    typedef typename op_table_t::tag_type tag_type;
    typedef boost::tuple<std::size_t, std::size_t, tag_type, typename Matrix::value_type> pretensor_value;
    typedef std::vector<pretensor_value> pretensor_t;

    op_table_ptr original_op_table = model.operators_table()->get_operator_table();
    op_table_ptr new_op_table( new op_table_t() );
    
    size_t length = lat.size();
    
    std::vector<tag_type> ident(lat.maximum_vertex_type()+1), fill(lat.maximum_vertex_type()+1);
    for (int type=0; type<ident.size(); ++type) {
        ident[type] = new_op_table->register_op( model.identity_matrix(type) );
        fill[type]  = new_op_table->register_op( model.filling_matrix(type)  );
    }
    
    MPO<Matrix, SymmGroup> mpo(length);
    std::vector<bool> used_p(length, false);
    
    for (int n=0; n<hamil_terms.size(); )
    {
        if (hamil_terms[n].size() != 2) throw std::runtime_error("hamiltonian terms have to be sorted with two-site bond terms before site terms.");
        int pos1 = hamil_terms[n].position(0);
        int pos2 = hamil_terms[n].position(1);
        
        int type1 = lat.get_prop<int>("type", pos1);
        int type2 = lat.get_prop<int>("type", pos2);
        
        exp_mpo_maker<Matrix, SymmGroup> maker(model.phys_dim(type1), model.phys_dim(type2), pos1, pos2, lat);
        
        maker.add_term(hamil_terms[n], model);
        
        int k = n+1;
        for (; k<hamil_terms.size() && hamil_terms[n].site_match(hamil_terms[k]); ++k)
            maker.add_term(hamil_terms[k], model);
        
        MPO<Matrix, SymmGroup> block_mpo = maker.exp_mpo(alpha, ident, fill, new_op_table);
        for (size_t p=0; p<pos2-pos1+1; ++p) {
            using std::swap;
            swap(mpo[pos1+p], block_mpo[p]);
            used_p[pos1+p] = true;
        }
        
        n = k;
    }
    
    
    // Filling missing identities
    for (std::size_t p=0; p<length; ++p) {
        if (!used_p[p]) {
            pretensor_t preident(1, pretensor_value(0, 0, ident[lat.get_prop<int>("type", p)], 1.) );
            MPOTensor<Matrix, SymmGroup> r(1, 1, preident, new_op_table);
            using std::swap;
            swap(mpo[p], r);
            used_p[p] = true;
        }
    }
    
    return mpo;
}

#endif
