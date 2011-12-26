/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuum/models_u1.hpp"

namespace app {
    
    template<class Matrix>
    struct cont_model_factory<Matrix, U1> {
        static typename model_traits<Matrix, U1>::model_ptr parse
        (Lattice const & lattice, BaseParameters & model)
        {
            std::string model_str = model.is_set("model") ? "model" : "MODEL";
            if (model.get<std::string>(model_str) == std::string("optical_lattice"))
                return typename model_traits<Matrix, U1>::model_ptr(
                            new OpticalLattice<Matrix>(lattice, model)
                       );
            else {
                throw std::runtime_error("Don't know this model!");
                return typename model_traits<Matrix, U1>::model_ptr();
            }
        }
    };

}
