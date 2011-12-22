/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuum/models_2u1.hpp"

namespace app {
    
    template<class Matrix>
    struct cont_model_factory<Matrix, TwoU1> {
        static typename model_traits<Matrix, TwoU1>::model_ptr parse
        (Lattice const & lattice, BaseParameters & model)
        {
            if (model.get<std::string>("MODEL") == std::string("fermi_optical_lattice"))
                return typename model_traits<Matrix, TwoU1>::model_ptr(
                            new FermiOpticalLattice<Matrix>(lattice, model)
                       );
            else {
                throw std::runtime_error("Don't know this model!");
                return typename model_traits<Matrix, TwoU1>::model_ptr();
            }
        }
    };

}
