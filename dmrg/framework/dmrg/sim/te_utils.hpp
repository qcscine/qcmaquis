
#ifndef APP_TE_UTILS_H
#define APP_TE_UTILS_H

#include "dmrg/models/hamiltonian.h"

#include "types/utils/matrix_vector_traits.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/multi_index.h"
#include "dmrg/mp_tensors/reshapes.h"
#include "dmrg/mp_tensors/generic_reshape.h"
#include "dmrg/mp_tensors/mpo_manip.h"

#include <vector>
#include <set>
#include <algorithm>

namespace app {
    
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
//        for (int n=0; n<ret.size(); ++n)
//        {
//            cout << "Hamiltonian #" << n << std::endl;
//            cout << ret[n];
//        }
        
        return ret;
    }
    
    // Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
    template <class Matrix, class SymmGroup>
    std::map<std::size_t, block_matrix<Matrix, SymmGroup> > make_exp_nn (Hamiltonian<Matrix, SymmGroup> const & H, typename Matrix::value_type const & alpha = 1)
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
    class term_exp_mpo {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef typename MultiIndex<SymmGroup>::index_id index_id;
        typedef typename MultiIndex<SymmGroup>::set_id set_id;

    public:
        term_exp_mpo(Index<SymmGroup>const & phys_, op_t const & ident_,
                     std::size_t pos1_, std::size_t pos2_,
                     bool s = true)
        : simple(s)
        , ident(ident_)
        , phys(phys_)
        , pos1(pos1_)
        , pos2(pos2_)
        {
            if (!simple) {
                std::vector<std::pair<index_id, bool> > left_vec, right_vec;
                for (size_t p=pos1; p<=pos2; ++p) {
                    index_id id1 = midx.insert_index(phys);
                    index_id id2 = midx.insert_index(phys);
                    
                    left_vec.push_back( std::make_pair(id1, true) );
                    right_vec.push_back( std::make_pair(id2, true) );
                    
                }
                op_set = midx.create_set(left_vec, right_vec);
            }
        }
        
        void add_term(Hamiltonian_Term<Matrix, SymmGroup> const & term)
        {
            op_t tmp;
            if (term.operators.size() == 2)
                kron_ops(term.operators[0].second, term.operators[1].second, term.fill_operator, tmp);
            else if (term.operators[0].first == pos1)
                kron_ops(term.operators[0].second, ident, ident, tmp);
            else if (term.operators[0].first == pos2)
                kron_ops(ident, term.operators[0].second, ident, tmp);
            else
                throw std::runtime_error("Operator k not matching any valid position.");
            bond_op += tmp;

        }
        
        inline MPO<Matrix, SymmGroup> get_exp_mpo(typename Matrix::value_type const & alpha = 1) const
        {
            if (simple)
                return get_simple_mpo(alpha);
            else
                return get_gen_mpo(alpha);
        }
        
    private:
        
        inline void kron_ops(op_t const &op1, op_t const & op2, op_t const & fill, op_t & bond) const
        {
            if (simple)
                op_kron(phys, op1, op2, bond);
            else
                op_kron_long(midx, op_set, op1, op2, fill, pos2-pos1, bond);
        }
        
        MPO<Matrix, SymmGroup> get_simple_mpo(typename Matrix::value_type const & alpha) const
        {
            op_t bond_exp;
            bond_exp = op_exp(phys*phys, bond_op, alpha);
            bond_exp = reshape_2site_op(phys, bond_exp);
            block_matrix<Matrix, SymmGroup> U, V, left, right;
            block_matrix<typename maquis::types::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Ssqrt;
            svd(bond_exp, U, V, S);
            Ssqrt = sqrt(S);
            gemm(U, Ssqrt, left);
            gemm(Ssqrt, V, right);
            
            // reshape and write back
            std::vector<block_matrix<Matrix, SymmGroup> > U_list = reshape_right_to_list(phys, left);
            std::vector<block_matrix<Matrix, SymmGroup> > V_list = reshape_left_to_list(phys, right);
            assert(U_list.size() == V_list.size());
            
            MPOTensor<Matrix, SymmGroup> left_tensor(1, U_list.size());
            MPOTensor<Matrix, SymmGroup> middle_tensor(U_list.size(), U_list.size());
            MPOTensor<Matrix, SymmGroup> right_tensor(U_list.size(), 1);
            
            for (std::size_t use_b=0; use_b<U_list.size(); ++use_b)
            {
                left_tensor(0, use_b) = U_list[use_b];
                middle_tensor(use_b, use_b) = ident;
                right_tensor(use_b, 0) = V_list[use_b];
            }
            
            MPO<Matrix, SymmGroup> mpo(pos2-pos1+1);
            mpo[0] = left_tensor;
            mpo[pos2-pos1] = right_tensor;
            for (std::size_t p=1; p<pos2-pos1; ++p)
            {
                mpo[p] = middle_tensor;
            }

            return mpo;
        }
        
        MPO<Matrix, SymmGroup> get_gen_mpo(typename Matrix::value_type const & alpha) const
        {
            Index<SymmGroup> op_basis = phys;
            for (size_t p=pos1+1; p<=pos2; ++p)
                    op_basis = op_basis * phys;
            
            op_t bond_exp = op_exp(op_basis, bond_op, alpha);
            MPO<Matrix, SymmGroup> block_mpo = block_to_mpo(phys, bond_exp, pos2-pos1+1);
            return block_mpo;
        }

        
        op_t ident;
        Index<SymmGroup> phys;
        std::size_t pos1, pos2;
        bool simple;
        op_t bond_op;
        
        MultiIndex<SymmGroup> midx;
        set_id op_set;
        
    };
    
    // Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
    template <class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_exp_mpo (std::size_t length,
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
            
            term_exp_mpo<Matrix, SymmGroup> t(H.get_phys(), H.get_identity(), pos1, pos2,
                                              !H[n].with_sign);
            
            t.add_term(H[n]);
            
            int k = n+1;
            for (; k<H.n_terms() && H[n].site_match(H[k]); ++k)
                t.add_term(H[k]);
            
            MPO<Matrix, SymmGroup> block_mpo = t.get_exp_mpo(alpha);
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
    
    
} // namespace

#endif
