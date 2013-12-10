/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/models_u1.hpp"
//#include "dmrg/models/coded/models_bela.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    static boost::shared_ptr<model_impl<Matrix, U1> > parse
    (Lattice const& lattice, BaseParameters & model)
    {
        typedef boost::shared_ptr<model_impl<Matrix, U1> > impl_ptr;
        if (model["MODEL"] == std::string("heisenberg"))
            return impl_ptr( new Heisenberg<Matrix>(lattice, model["Jxy"], model["Jz"]) );
        else if (model["MODEL"] == std::string("HCB"))
            return impl_ptr( new HCB<Matrix>(lattice) );
        else if (model["MODEL"] == std::string("boson Hubbard"))
            return impl_ptr( new BoseHubbard<Matrix>(lattice, model) );
//        else if (model["MODEL"] == std::string("fermion Hubbard"))
//            return impl_ptr( new FermiHubbardU1<Matrix>(lattice, model) );
        else if (model["MODEL"] == std::string("FreeFermions"))
            return impl_ptr( new FreeFermions<Matrix>(lattice, model["t"]) );
//        else if (model["MODEL"] == std::string("bela_chiral"))
//            return impl_ptr( new Chiral<Matrix>(lattice, model) );
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
