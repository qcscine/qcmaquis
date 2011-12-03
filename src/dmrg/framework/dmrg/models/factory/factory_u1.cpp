/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "factory_impl.h"

namespace app {
	typedef U1 grp;

#ifdef USE_MTM
    impl_init_model(mtmatrix1, grp)
    impl_init_model(cmtmatrix1, grp)
#else
    impl_init_model(matrix1, grp)
    impl_init_model(cmatrix1, grp)
#endif


    template<class Matrix>
    struct hamil_factory<Matrix, U1> {
        static Hamiltonian<Matrix, U1> parse (BaseParameters & model, Lattice const & lattice)
        {
            if (model.get<std::string>("MODEL") == std::string("heisenberg"))
                return Heisenberg<Matrix>(lattice, model.get<double>("Jxy"), model.get<double>("Jz"));
            else if (model.get<std::string>("MODEL") == std::string("HCB"))
                return HCB<Matrix>(lattice);
            else if (model.get<std::string>("MODEL") == std::string("FreeFermions"))
                return FreeFermions<Matrix>(lattice, model.get<double>("t"));
            else {
                throw std::runtime_error("Don't know this model!");
                return Hamiltonian<Matrix, U1>();
            }
        }
    };
    
    template <>
    U1::charge init_qn<U1> (BaseParameters & model)
    {
        return model.get<int>("u1_total_charge");
    }
}
