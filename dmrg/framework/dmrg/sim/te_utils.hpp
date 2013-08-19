
#ifndef APP_TE_UTILS_H
#define APP_TE_UTILS_H

#include "dmrg/models/hamiltonian.h"

#include "utils/traits.hpp"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/generic_reshape.h"
#include "dmrg/mp_tensors/mpo_manip.h"
#include "dmrg/mp_tensors/compression.h"


#include <vector>
#include <set>
#include <algorithm>

// Return: Bond terms are allowed to be in the same Hamiltonian object if they do not overlap
//         Site terms are splitted among all Hamiltonian objects using that site
template <class Matrix, class SymmGroup>
std::vector<Hamiltonian<Matrix, SymmGroup> > separate_overlaps (Hamiltonian<Matrix, SymmGroup> const & H)
{
    typedef std::map<std::size_t, std::set<std::size_t> > pos_where_t;
    typedef std::map<std::size_t, std::vector<typename Hamiltonian<Matrix, SymmGroup>::hamterm_t> > pos_terms_t;
    
    std::vector<Hamiltonian<Matrix, SymmGroup> > ret;
    pos_where_t pos_where;
    pos_terms_t pos_terms;
    
    for (int i=0; i<H.n_terms(); ++i)
    {            
        if (H[i].operators.size() == 1) {
            pos_terms[H[i].operators[0].first].push_back(H[i]);
            continue;
        }
        
        Hamiltonian_Term<Matrix, SymmGroup> term = H[i];
        term.canonical_order(); // TODO: check and fix for fermions!!
        bool used = false;
        for (int n=0; n<ret.size() && !used; ++n)
        {
            bool overlap = false;
            for (int j=0; j<ret[n].n_terms() && !overlap; ++j)
            {
                if ( ret[n][j].site_match(term) ) break;
                overlap = ret[n][j].overlap(term);
            }
            
            if (!overlap) {
            	ret[n].add_term(term);
            	for (int p=0; p<term.operators.size(); ++p)
            		pos_where[term.operators[p].first].insert(n);
                used = true;
            }
        }
        
        if (!used) {
            Hamiltonian<Matrix, SymmGroup> tmp;
            tmp.set_identity( H.get_identity() );
            tmp.set_phys( H.get_phys() );
            tmp.add_term(term);
            ret.push_back(tmp);
            
        	for (int p=0; p<term.operators.size(); ++p)
        		pos_where[term.operators[p].first].insert(ret.size()-1);
        }
    }
    
    // Adding site terms to all Hamiltonians acting on site i
    for (typename pos_terms_t::const_iterator it = pos_terms.begin();
         it != pos_terms.end();
         ++it)
    {
        double coeff = 1. / pos_where[it->first].size();
        
        for (typename std::set<std::size_t>::const_iterator it2 = pos_where[it->first].begin();
             it2 != pos_where[it->first].end();
             ++it2)
            for (int k=0; k<it->second.size(); ++k)
            {
            	typename Hamiltonian<Matrix, SymmGroup>::hamterm_t tmp_term = it->second[k];
            	tmp_term.operators[0].second *= coeff;
                ret[*it2].add_term(tmp_term);
            }
    }
    
    // Sorting individual Hamiltonians
    for (int n=0; n<ret.size(); ++n)
        std::sort(ret[n].begin(), ret[n].end());
    
    // Output
    //    for (int n=0; n<ret.size(); ++n)
    //    {
    //        maquis::cout << "Hamiltonian #" << n << std::endl;
    //        maquis::cout << ret[n];
    //    }
    
    return ret;
}

// Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
template <class Matrix, class SymmGroup>                                                                                                // (todo: 30.04.12 / Matthias scalar/value types discussion)
std::map<std::size_t, block_matrix<Matrix, SymmGroup> > make_exp_nn (Hamiltonian<Matrix, SymmGroup> const & H, typename Matrix::value_type const & alpha = 1.) // type of the time step // template
{
    typedef Hamiltonian<Matrix, SymmGroup> ham;
    
    std::map<std::size_t, block_matrix<Matrix, SymmGroup> > map_exp;
    
    for (int n=0; n<H.n_terms(); )
    {
        assert(H[n].operators.size() == 2);
        int pos1 = H[n].operators[0].first;
        int pos2 = H[n].operators[1].first;
        assert(std::abs(pos1-pos2) == 1);
        
        typename ham::op_t bond_op;
        op_kron(H.get_phys(), H[n].operators[0].second, H[n].operators[1].second, bond_op);
        
        int k = n+1;
        for (; k<H.n_terms(); ++k)
        {
            Hamiltonian_Term<Matrix, SymmGroup> term = H[k];
            if (! H[n].site_match(term))
                break;
            typename ham::op_t tmp;
            if (term.operators.size() == 2)
                op_kron(H.get_phys(), term.operators[0].second, term.operators[1].second, tmp);
            else if (term.operators[0].first == pos1)
                op_kron(H.get_phys(), term.operators[0].second, H.get_identity(), tmp);
            else if (term.operators[0].first == pos2)
                op_kron(H.get_phys(), H.get_identity(), term.operators[0].second, tmp);
            else
                throw std::runtime_error("Operator k not matching any valid position.");
            bond_op += tmp;
        }
        
        bond_op = op_exp(H.get_phys()*H.get_phys(), bond_op, alpha);
        
        map_exp[pos1] = bond_op;
        
        n = k;
    }
            
    return map_exp;
}


template <class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> term2block(typename Hamiltonian<Matrix, SymmGroup>::hamterm_t const & term, std::size_t pos1,
                                           Index<SymmGroup> const & phys_i, block_matrix<Matrix, SymmGroup> const & ident)
{
#ifndef NDEBUG
    if (term.operators.size() == 2)
        assert(std::abs( term.operators[0].first-term.operators[1].first ) == 1);
#endif
    
    block_matrix<Matrix, SymmGroup> bond_op;
    if (term.operators.size() == 2)
        op_kron(phys_i, term.operators[0].second, term.operators[1].second, bond_op);
    else if (term.operators[0].first == pos1)
        op_kron(phys_i, term.operators[0].second, ident, bond_op);
    else
        op_kron(phys_i, ident, term.operators[0].second, bond_op);
//    else
//        throw std::runtime_error("Operator k not matching any valid position.");
    return bond_op;
}

template <class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > hamil_to_blocks(Hamiltonian<Matrix, SymmGroup> const & H, std::size_t L)
{
    std::vector<block_matrix<Matrix, SymmGroup> > ret_blocks(L-1);
    
    for (int i=0; i<H.n_terms(); ++i)
    {
        Hamiltonian_Term<Matrix, SymmGroup> term = H[i];
        term.canonical_order(); // TODO: check and fix for fermions!!
        std::size_t pos1 = term.operators[0].first;
        if (term.operators.size() == 1) {
            if (pos1 != 0 && pos1 != L-1)
                term.operators[0].second /= 2.;
            if (pos1 < L-1)
                ret_blocks[pos1] += term2block(term, pos1, H.get_phys(), H.get_identity());
            if (pos1 > 0)
                ret_blocks[pos1-1] += term2block(term, pos1-1, H.get_phys(), H.get_identity());
        } else if (term.operators.size() == 2) {
            ret_blocks[pos1] += term2block(term, pos1, H.get_phys(), H.get_identity());
        }
    }
    
    return ret_blocks;
}

template <class Matrix, class SymmGroup>
class exp_mpo_maker {
    typedef block_matrix<Matrix, SymmGroup> op_t;

public:
    exp_mpo_maker(Index<SymmGroup> const& phys_, op_t const& ident_,
                 std::size_t pos1_, std::size_t pos2_)
    : ident(ident_)
    , phys(phys_)
    , pos1(pos1_), pos2(pos2_), L(pos2-pos1+1)
    , n_boso(0), n_ferm(0)
    { }
    
    void add_term(Hamiltonian_Term<Matrix, SymmGroup> const & term, typename Matrix::value_type const & alpha = 1.)
    {
        if (term.operators.size() > 2)
            throw std::runtime_error("time evolution requires at max bond term.");
        
        /// kron product of operators
        op_t bond_op;
        if (term.operators.size() == 2)
            op_kron(phys, term.operators[0].second, term.operators[1].second, bond_op);
        else if (term.operators[0].first == pos1)
            op_kron(phys, term.operators[0].second, ident, bond_op);
        else if (term.operators[0].first == pos2)
            op_kron(phys, ident, term.operators[0].second, bond_op);
        else
            throw std::runtime_error("Operator k not matching any valid position.");
        
        if (term.with_sign) {
            fermionic_bond += bond_op;
            n_ferm += 1;
            fill = term.fill_operator;
        } else {
            bosonic_bond += bond_op;
            n_boso += 1;
        }
    }
    
    MPO<Matrix, SymmGroup> exp_mpo(typename Matrix::value_type const & alpha = 1.) const
    {
        std::size_t maximum_b = 0;
        MPO<Matrix, SymmGroup> mpo(L);
        
        if (n_boso > 0) {
            std::vector<op_t> left_ops, right_ops;
            exp_and_split(bosonic_bond, alpha, left_ops, right_ops);
            
            maximum_b = add_to_mpo(mpo, maximum_b, left_ops, right_ops, 1., ident);
        }
        
        if (n_ferm > 0) {
            /// exp(alpha * op1 \otimes op2)
            {
                std::vector<op_t> left_ops, right_ops;
                exp_and_split(fermionic_bond, alpha, left_ops, right_ops);
                
                maximum_b = add_to_mpo(mpo, maximum_b, left_ops, right_ops, .5, ident);
                maximum_b = add_to_mpo(mpo, maximum_b, left_ops, right_ops, .5, fill);
            }
            
            /// exp(-alpha * op1 \otimes op2)
            {
                std::vector<op_t> left_ops, right_ops;
                exp_and_split(fermionic_bond, -alpha, left_ops, right_ops);
                
                maximum_b = add_to_mpo(mpo, maximum_b, left_ops, right_ops, .5, ident);
                maximum_b = add_to_mpo(mpo, maximum_b, left_ops, right_ops, -.5, fill);
            }
        }
        return mpo;
    }
    
private:
    
    void exp_and_split(op_t const& bond_op, typename Matrix::value_type const & alpha,
                       std::vector<op_t> & left_ops, std::vector<op_t> &  right_ops) const
    {
        op_t bond_exp;
        bond_exp = op_exp_hermitian(phys*phys, bond_op, alpha);
        bond_exp = reshape_2site_op(phys, bond_exp);
        typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type diag_matrix;
        block_matrix<Matrix, SymmGroup> U, V, left, right;
        block_matrix<diag_matrix, SymmGroup> S, Ssqrt;
        svd(bond_exp, U, V, S);
        Ssqrt = sqrt(S);
        gemm(U, Ssqrt, left);
        gemm(Ssqrt, V, right);
        
        for(std::size_t k = 0; k < S.n_blocks(); ++k){
            std::vector<typename diag_matrix::value_type> sk = maquis::bindings::matrix_cast< std::vector<typename diag_matrix::value_type> >(S[k]);
            int keep = std::find_if(sk.begin(), sk.end(), boost::lambda::_1 < 1e-10)-sk.begin();
            
            left.resize_block(left.left_basis()[k].first, left.right_basis()[k].first,
                              left.left_basis()[k].second, keep);
            right.resize_block(right.left_basis()[k].first, right.right_basis()[k].first,
                               keep, right.right_basis()[k].second);
        }
        
        left_ops  = reshape_right_to_list(phys, left);
        right_ops = reshape_left_to_list(phys, right);
        assert(left_ops.size() == right_ops.size());
    }
    
    std::size_t add_to_mpo(MPO<Matrix, SymmGroup> & mpo, std::size_t maximum_b, 
                           std::vector<op_t> const& left_ops, std::vector<op_t> const& right_ops,
                           double s, op_t const& fill) const
    {
        for (std::size_t i=0; i<left_ops.size(); ++i)
        {
            std::size_t b = (maximum_b++);
            mpo[0](0, b) = s*left_ops[i];
            mpo[L-1](b, 0) = right_ops[i];
            for (std::size_t p=1; p<L-1; ++p)
                mpo[p](b, b) = fill;
        }
        return maximum_b;
    }
    
    op_t ident, fill;
    Index<SymmGroup> phys;
    std::size_t pos1, pos2, L;
    std::size_t n_boso, n_ferm;
    op_t bosonic_bond, fermionic_bond;
};

// Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
template <class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_exp_mpo(std::size_t length,
                                    Hamiltonian<Matrix, SymmGroup> const & H,
                                    typename Matrix::value_type const & alpha = 1)
{
    typedef Hamiltonian<Matrix, SymmGroup> ham;
    
    MPO<Matrix, SymmGroup> mpo(length);
    std::vector<bool> used_p(length, false);
    
    for (int n=0; n<H.n_terms(); )
    {
        assert(H[n].operators.size() == 2);
        int pos1 = H[n].operators[0].first;
        int pos2 = H[n].operators[1].first;
        
        exp_mpo_maker<Matrix, SymmGroup> maker(H.get_phys(), H.get_identity(), pos1, pos2);
        
        maker.add_term(H[n]);
        
        int k = n+1;
        for (; k<H.n_terms() && H[n].site_match(H[k]); ++k)
            maker.add_term(H[k]);
        
        MPO<Matrix, SymmGroup> const& block_mpo = maker.exp_mpo(alpha);
        for (size_t p=0; p<pos2-pos1+1; ++p) {
            mpo[pos1+p] = block_mpo[p];
            used_p[pos1+p] = true;
        }
        
        n = k;
    }
    
    // Filling missing identities
    for (std::size_t p=0; p<length; ++p)
        if (!used_p[p]) {
            MPOTensor<Matrix, SymmGroup> r(1, 1);
            r(0, 0) = H.get_identity();
            mpo[p] = r;
            used_p[p] = true;
        }
    
    return mpo;
}

#endif
