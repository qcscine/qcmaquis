/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/models_u1.hpp"
#include "dmrg/models/coded/models_bela.hpp"

template<class Matrix>
struct model_factory<Matrix, U1> {
    static typename model_traits<Matrix, U1>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        if (model["MODEL"] == std::string("heisenberg"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new Heisenberg<Matrix>(lattice, model["Jxy"], model["Jz"])
                   );
        else if (model["MODEL"] == std::string("HCB"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new HCB<Matrix>(lattice)
                   );
        else if (model["MODEL"] == std::string("boson Hubbard"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new BoseHubbard<Matrix>(lattice, model["Nmax"], model["t"], model["U"], model["V"])
                                                                );
        else if (model["MODEL"] == std::string("fermion Hubbard"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new FermiHubbardU1<Matrix>(lattice, model)
                                                                );
        else if (model["MODEL"] == std::string("FreeFermions"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new FreeFermions<Matrix>(lattice, model["t"])
                   );
        else if (model["MODEL"] == std::string("bela_chiral"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new Chiral<Matrix>(lattice, model)
                   );
        else {
            throw std::runtime_error("Don't know this model!");
            return typename model_traits<Matrix, U1>::model_ptr();
        }
    }
};
