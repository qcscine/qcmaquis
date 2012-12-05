/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/models_none.hpp"

template<class Matrix>
struct model_factory<Matrix, TrivialGroup> {
    static typename model_traits<Matrix, TrivialGroup>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        if (model.get<std::string>("MODEL") == std::string("boson Hubbard"))
            return typename model_traits<Matrix, TrivialGroup>::model_ptr(
                                                                          new BoseHubbardNone<Matrix>(lattice, model.get<int>("Nmax"), model.get<double>("t"), model.get<double>("U"), model.get<double>("V"))
                                                                );
        else {
            throw std::runtime_error("Don't know this model with NONE symmetry group!");
            return typename model_traits<Matrix, TrivialGroup>::model_ptr();
        }
    }
};
