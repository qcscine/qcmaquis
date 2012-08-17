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
        if (model.get<std::string>("MODEL") == std::string("heisenberg"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new Heisenberg<Matrix>(lattice, model.get<double>("Jxy"), model.get<double>("Jz"))
                   );
        else if (model.get<std::string>("MODEL") == std::string("HCB"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new HCB<Matrix>(lattice)
                   );
        else if (model.get<std::string>("MODEL") == std::string("boson Hubbard"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new BoseHubbard<Matrix>(lattice, model.get<int>("Nmax"), model.get<double>("t"), model.get<double>("U"), model.get<double>("V"))
                                                                );
        else if (model.get<std::string>("MODEL") == std::string("fermion Hubbard"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new FermiHubbardU1<Matrix>(lattice, model)
                                                                );
        else if (model.get<std::string>("MODEL") == std::string("FreeFermions"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new FreeFermions<Matrix>(lattice, model.get<double>("t"))
                   );
        else if (model.get<std::string>("MODEL") == std::string("bela_chiral"))
            return typename model_traits<Matrix, U1>::model_ptr(
                        new Chiral<Matrix>(lattice, model)
                   );
        else {
            throw std::runtime_error("Don't know this model!");
            return typename model_traits<Matrix, U1>::model_ptr();
        }
    }
};
