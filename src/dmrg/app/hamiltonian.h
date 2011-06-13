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

#include "generate_mpo.hpp"

namespace app {
    
    enum TermsType {All_terms, Site_terms, Odd_terms, Even_terms, ExpSite_terms, ExpOdd_terms, ExpEven_terms};
    
    template<class Matrix, class SymmGroup>
    struct Hamiltonian_Term : public generate_mpo::Operator_Term<Matrix, SymmGroup>
    {
    	typedef typename generate_mpo::Operator_Term<Matrix, SymmGroup>::op_t op_t;
    };
    
    template<class Matrix, class SymmGroup>
    class Hamiltonian
    {
    public:
        typedef Hamiltonian_Term<Matrix, SymmGroup> hamterm_t;
        typedef typename hamterm_t::op_t op_t;
        typedef typename std::vector<hamterm_t>::const_iterator const_iterator;
        typedef typename std::vector<hamterm_t>::iterator iterator;
        
        Hamiltonian () {}
        Hamiltonian (Index<SymmGroup> const & phys_,
                     op_t const & ident_,
                     std::vector<hamterm_t> const & terms_)
        : phys(phys_)
        , ident(ident_)
        , terms(terms_)
        {}
        
        virtual int n_terms (TermsType what=All_terms) const { return terms.size(); }
        virtual hamterm_t const & operator[] (int i) const { return terms[i]; }
        virtual void add_term (hamterm_t const & term) { terms.push_back(term); }
        
        virtual Index<SymmGroup> get_phys () const { return phys; }
        virtual void set_phys (Index<SymmGroup> const & phys_) { phys = phys_; }
        virtual op_t const & get_identity () const { return ident; }
        virtual void set_identity (op_t const & ident_) { ident = ident_; }
        
        const_iterator begin () const { return terms.begin(); }
        const_iterator end () const { return terms.end(); }
        
        iterator begin () { return terms.begin(); }
        iterator end () { return terms.end(); }

    protected:
        std::vector<hamterm_t> terms;
        Index<SymmGroup> phys;
        op_t ident;
    };
}

// call me to do something!

namespace app {
    
    
    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_mpo(std::size_t L, Hamiltonian<Matrix, SymmGroup> const & H, TermsType what=All_terms)
    {
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(L, H.get_identity());
        
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
