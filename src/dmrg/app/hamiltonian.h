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

#include "lattice.h"

namespace app {
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
        virtual int n_terms() const = 0;
        virtual Hamiltonian_Term<Matrix, SymmGroup> operator[](int) const = 0;
        virtual Index<SymmGroup> get_phys() const = 0;
        virtual typename Hamiltonian_Term<Matrix, SymmGroup>::op_t get_identity() const = 0;
    };
}

#include "hamiltonian_detail.hpp"

// call me to do something!

namespace app {
    template<class Matrix, class SymmGroup>
    MPO<Matrix, SymmGroup> make_mpo(std::size_t L, Hamiltonian<Matrix, SymmGroup> const & H)
    {
        hamiltonian_detail::MPOMaker<Matrix, SymmGroup> mpom(L, H.get_identity());
        
        for (int i = 0; i < H.n_terms(); ++i)
        {
            mpom.add_term(H[i]);
        }
        
        return mpom.create_mpo();
    }
}

#endif
