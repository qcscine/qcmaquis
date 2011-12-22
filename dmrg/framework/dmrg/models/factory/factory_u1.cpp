/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

#include "dmrg/models/coded/factory_u1.hpp"
#include "dmrg/models/continuum/factory_u1.hpp"

namespace app {
	typedef U1 grp;

#ifdef USE_MTM
    impl_init_model(mtmatrix1, grp)
    impl_init_model(cmtmatrix1, grp)
#else
    impl_init_model(matrix1, grp)
    impl_init_model(cmatrix1, grp)
#endif

    
    template <>
    U1::charge init_qn<U1> (BaseParameters & model)
    {
        return model.get<int>("u1_total_charge");
    }
}
