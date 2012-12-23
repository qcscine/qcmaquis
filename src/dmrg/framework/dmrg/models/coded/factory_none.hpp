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
struct model_factory<Matrix, TrivialGroup> {
    static typename model_traits<Matrix, TrivialGroup>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        if (model.get<std::string>("MODEL") == std::string("boson Hubbard"))
            return typename model_traits<Matrix, TrivialGroup>::model_ptr(
                                                                          new BoseHubbardNone<Matrix>(lattice, model.get<int>("Nmax"), model.get<double>("t"), model.get<double>("U"), model.get<double>("V"))
                                                                );
        else if (model.get<std::string>("MODEL") == std::string("super boson Hubbard"))
            return typename model_traits<Matrix, TrivialGroup>::model_ptr(
                                                                          new SuperBoseHubbardNone<Matrix>(lattice, model.get<int>("Nmax"), model.get<double>("t"), model.get<double>("U"), model.get<double>("V"), model.get<double>("Delta"), model.get<double>("Gamma1a"), model.get<double>("Gamma1b"), model.get<double>("Gamma2"))
                                                                );
        else {
            throw std::runtime_error("Don't know this model with NONE symmetry group!");
            return typename model_traits<Matrix, TrivialGroup>::model_ptr();
        }
    }
};
