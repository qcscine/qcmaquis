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
struct model_factory<Matrix, Ztwo> {
    static boost::shared_ptr<model_impl<Matrix, SymmGroup> > parse
    (Lattice const & lattice, BaseParameters & model)
    { }
};
