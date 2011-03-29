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

#include "lattice.h"

namespace app {
    template<class Matrix, class SymmGroup>
    struct Hamiltonian_Term
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        
        std::vector<std::pair<typename Lattice::pos_t, op_t> > operators;
        op_t fill_operator;
    }
    
    // implement me for your purpose!
    template<class Matrix, class SymmGroup>
    class Hamiltonian
    {
    public:
        virtual int n_terms() const = 0;
        virtual Hamiltonian_Term<Matrix, SymmGroup> operator[](int) const = 0;
    };
}

#include "hamiltonian_detail.hpp"

// call me to do something!

namespace app {
    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_mpo(Lattice const & adj,
                                    Hamiltonian<Matrix, SymmGroup> const & H)
    {
        MPOMaker<Matrix, SymmGroup> mpom;
        
        for (int i = 0; i < H.n_terms(); ++i)
        {
            mpom.add_term(H[i]);
        }
        
        return mpom.create_mpo();
    }
}

#endif
