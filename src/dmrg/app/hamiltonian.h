/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef HAMILTONIAN_H
#define HAMILTONIAN_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/block_matrix_algorithms.h"
#include "block_matrix/symmetry.h"

#include <vector>
#include <iostream>

#include "lattice.h"

namespace app {
    
    enum TermsType {All_terms, Site_terms, Odd_terms, Even_terms, ExpSite_terms, ExpOdd_terms, ExpEven_terms};

    template<class Matrix, class SymmGroup>
    struct Hamiltonian_Term
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        
        std::vector<std::pair<typename Lattice::pos_t, op_t> > operators;
        op_t fill_operator;
    };
    
    // implement me for your purpose!
    template<class Matrix, class SymmGroup>
    class Hamiltonian
    {
    public:
        virtual int n_terms(TermsType what=All_terms) const = 0;
        virtual Hamiltonian_Term<Matrix, SymmGroup> operator[](int) const = 0;

        /*
        virtual int n_site_terms() const = 0;
        virtual Hamiltonian_Term<Matrix, SymmGroup> site_term(int) const = 0;
        virtual int n_odd_bonds() const = 0; // 1->2, 3->4, ...
        virtual Hamiltonian_Term<Matrix, SymmGroup> odd_bond(int) const = 0;
        virtual int n_even_bonds() const = 0; // 2->3, 4->5, ...
        virtual Hamiltonian_Term<Matrix, SymmGroup> even_bond(int) const = 0;
        
        // TODO: implement default behaviour = exp+svd to decompose
        virtual Hamiltonian_Term<Matrix, SymmGroup> exp_site_term(int i) const
        {
            return Hamitonian_Term<Matrix>();
        }
        virtual Hamiltonian_Term<Matrix, SymmGroup> exp_odd_bond(int i) const
        {
            return Hamitonian_Term<Matrix>();
        }
        virtual Hamiltonian_Term<Matrix, SymmGroup> exp_even_bond(int i) const
        {
            return Hamitonian_Term<Matrix>();
        }
        */

        virtual Index<SymmGroup> get_phys() const = 0;
        virtual typename Hamiltonian_Term<Matrix, SymmGroup>::op_t get_identity() const = 0;
    };
}

#include "hamiltonian_detail.hpp"

// call me to do something!

namespace app {
    
    
    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_mpo(std::size_t L, Hamiltonian<Matrix, SymmGroup> const & H, TermsType what=All_terms)
    {
        hamiltonian_detail::MPOMaker<Matrix, SymmGroup> mpom(L, H.get_identity());
        
        for (int i = 0; i < H.n_terms(what); ++i)
        {
            mpom.add_term(H[i]);
        }
        
        return mpom.create_mpo();
    }
}


template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, app::Hamiltonian_Term<Matrix, SymmGroup> const & h)
{
    os << " - Fill operator:" << std::endl << h.fill_operator;
    for (int i=0; i<h.operators.size(); ++i)
        os << " - Operator at " << h.operators[i].first << ":" << std::endl << h.operators[i].second;
    return os;
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, app::Hamiltonian<Matrix, SymmGroup> const & H)
{
    for (int i=0; i<H.n_terms(); ++i)
        os << "* TERM(" << i << ") *" << std::endl << H[i];
    return os;
}


#endif
