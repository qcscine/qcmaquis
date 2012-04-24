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
        throw std::runtime_error("No models with TrivialGroup defined in the factory!");
        return typename model_traits<Matrix, TrivialGroup>::model_ptr();
    }
};
