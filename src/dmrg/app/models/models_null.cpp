/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "models_impl.h"

namespace app {
	typedef NullGroup grp;

#ifdef USE_MTM
    impl_init_model(mtmatrix1, grp)
#else
    impl_init_model(matrix1, grp)
    // impl_init_model(matrix2, grp)
    impl_init_model(cmatrix1, grp)
    // impl_init_model(cmatrix2, grp)
#endif


    template<class Matrix>
    struct hamil_factory<Matrix, NullGroup> {
        static Hamiltonian<Matrix, NullGroup> parse (BaseParameters & model, Lattice const & lattice)
        {
            throw std::runtime_error("No models with NullGroup defined in the factory!");
            return Hamiltonian<Matrix, NullGroup>();
        }
    };
    
    template <>
    NullGroup::charge init_qn<NullGroup> (BaseParameters & model)
    {
        return NullGroup::SingletCharge;
    }

}
