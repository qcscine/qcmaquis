/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "models_impl.h"

namespace app {
	typedef TwoU1 grp;

#ifdef USE_MTM
    impl_init_model(mtmatrix1, grp)
#else
    impl_init_model(matrix1, grp)
    // impl_init_model(matrix2, grp)
    impl_init_model(cmatrix1, grp)
    // impl_init_model(cmatrix2, grp)
#endif


    template<class Matrix>
    struct hamil_factory<Matrix, TwoU1> {
        static Hamiltonian<Matrix, TwoU1> parse (BaseParameters & model, Lattice const & lattice)
        {
            if (model.get<std::string>("model") == std::string("fermi_hubbard"))
                return TwoU1_FermiHubbard<Matrix>(lattice, model)();
            else {
                throw std::runtime_error("Don't know this model!");
                return Hamiltonian<Matrix, TwoU1>();
            }
        }
    };
    
    template <>
    TwoU1::charge init_qn<TwoU1> (BaseParameters & model)
    {
        TwoU1::charge initc;
        initc[0] = model.get<int>("u1_total_charge1");
        initc[1] = model.get<int>("u1_total_charge2");
        return initc;
    }


}
