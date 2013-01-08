/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuum/models_none.hpp"

template<class Matrix>
struct cont_model_factory<Matrix, Ztwo> {
    static typename model_traits<Matrix, Ztwo>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        throw std::runtime_error("Don't know this model!");
        return typename model_traits<Matrix, Ztwo>::model_ptr();
    }
};
