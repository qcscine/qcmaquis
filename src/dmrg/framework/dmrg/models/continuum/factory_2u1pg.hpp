/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuum/models_2u1.hpp"
#include "dmrg/models/continuum/super_models_2u1.hpp"

template<class Matrix>
struct cont_model_factory<Matrix, TwoU1PG> {
    static typename model_traits<Matrix, TwoU1PG>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        throw std::runtime_error("No Continuum models for TwoU1PG available!\n");
        return typename model_traits<Matrix, TwoU1PG>::model_ptr();
    }
};
