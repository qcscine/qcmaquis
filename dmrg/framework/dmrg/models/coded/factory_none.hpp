/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2012      by Jan Gukelberger <gukelberger@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/models_none.hpp"
#include "dmrg/models/coded/super_models_none.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, TrivialGroup> {
    static boost::shared_ptr<model_impl<Matrix, TrivialGroup> > parse
    (Lattice const& lattice, BaseParameters & model)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TrivialGroup> > impl_ptr;
        if (model["MODEL"] == std::string("boson Hubbard"))
            return impl_ptr( new BoseHubbardNone<Matrix>(lattice, model) );
        else if (model["MODEL"] == std::string("super boson Hubbard"))
            return impl_ptr( new SuperBoseHubbardNone<Matrix>(lattice, model) );
        else {
            throw std::runtime_error("Don't know this model with NONE symmetry group!");
            return impl_ptr();
        }
    }
};
