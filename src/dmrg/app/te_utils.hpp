
#ifndef APP_TE_UTILS_H
#define APP_TE_UTILS_H

#include "hamiltonian.h"

#include "dense_matrix/matrix_vector_traits.h"
#include "mp_tensors/reshapes.h"

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
            
            bool used = false;
            for (int n=0; n<ret.size() && !used; ++n)
            {
                bool overlap = false;
                for (int j=0; j<ret[n].n_terms() && !overlap; ++j)
                {
                    if ( ret[n][j].site_match(H[i]) ) break;
                    overlap = ret[n][j].overlap(H[i]);
                }
                
                if (!overlap) {
                	ret[n].add_term(H[i]);
                	for (int p=0; p<H[i].operators.size(); ++p)
                		pos_where[H[i].operators[p].first].insert(n);
                    used = true;
                }
            }
            
            if (!used) {
                Hamiltonian<Matrix, SymmGroup> tmp;
                tmp.set_identity( H.get_identity() );
                tmp.set_phys( H.get_phys() );
                tmp.add_term(H[i]);
                ret.push_back(tmp);
                
            	for (int p=0; p<H[i].operators.size(); ++p)
            		pos_where[H[i].operators[p].first].insert(ret.size()-1);
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
        for (int n=0; n<ret.size(); ++n)
        {
            std::cout << "Hamiltonian #" << n << std::endl;
            std::cout << ret[n];
        }
        
        return ret;
    }
    
    // Precondition: Hamiltonian has to be sorted with bond terms coming before site terms (default behaviour of Operator_Term::operator<())
    template <class Matrix, class SymmGroup>
    Hamiltonian<Matrix, SymmGroup> exp (Hamiltonian<Matrix, SymmGroup> const & H, typename Matrix::value_type const & alpha = 1)
    {
        typedef Hamiltonian<Matrix, SymmGroup> ham;
        
        
        for (int n=0; n<H.n_terms(); )
        {
            std::cout << "new group starting at n=" << n << std::endl;
            assert(H[n].operators.size() == 2);
            int pos1 = H[n].operators[0].first;
            int pos2 = H[n].operators[1].first;
            typename ham::op_t bond_op;
            op_kron(H.get_phys(), H[n].operators[0].second, H[n].operators[1].second, bond_op);
            
            int k = n+1;
            for (; k<H.n_terms() && H[n].site_match(H[k]); ++k)
            {
                std::cout << "using k=" << k << std::endl;
                typename ham::op_t tmp;
                if (H[k].operators.size() == 2)
                    op_kron(H.get_phys(), H[k].operators[0].second, H[k].operators[1].second, tmp);
                else if (H[k].operators[0].first == pos1)
                    op_kron(H.get_phys(), H[k].operators[0].second, H.get_identity(), tmp);
                else if (H[k].operators[0].first == pos2)
                    op_kron(H.get_phys(), H.get_identity(), H[k].operators[0].second, tmp);
                else
                    std::runtime_error("Operator k not matching any valid position.");
                bond_op += tmp;
            }
            std::cout << "group finishing with k=" << k << std::endl;
            
            bond_op = op_exp(H.get_phys()*H.get_phys(), bond_op, alpha);
            bond_op = reshape_2site_op(H.get_phys(), bond_op);
            block_matrix<Matrix, SymmGroup> U, V, left, right;
            block_matrix<typename blas::associated_diagonal_matrix<Matrix>::type, SymmGroup> S, Ssqrt;
            svd(bond_op, U, V, S);
            Ssqrt = sqrt(S);
            gemm(U, Ssqrt, left);
            gemm(Ssqrt, V, right);
            
            // reshape and write back
            std::vector<block_matrix<Matrix, SymmGroup> > U_list = reshape_right_to_list(H.get_phys(), left);
            std::vector<block_matrix<Matrix, SymmGroup> > V_list = reshape_left_to_list(H.get_phys(), right);
            
            
            n = k;
        }
        
        return H;
    }
    
} // namespace

#endif
