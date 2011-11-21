/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

namespace app {
	typedef TrivialGroup grp;

#ifdef USE_MTM
    impl_init_model(mtmatrix1, grp)
#else
    impl_init_model(matrix1, grp)
    impl_init_model(cmatrix1, grp)
#endif


    template<class Matrix>
    struct hamil_factory<Matrix, TrivialGroup> {
        static Hamiltonian<Matrix, TrivialGroup> parse (BaseParameters & model, Lattice const & lattice)
        {
            throw std::runtime_error("No models with TrivialGroup defined in the factory!");
            return Hamiltonian<Matrix, TrivialGroup>();
        }
    };
    
    template <>
    TrivialGroup::charge init_qn<TrivialGroup> (BaseParameters & model)
    {
        return TrivialGroup::IdentityCharge;
    }

}
