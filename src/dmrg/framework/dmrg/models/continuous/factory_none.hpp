/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/continuous/models_none.hpp"

namespace app {
    
    template<class Matrix>
    struct cont_model_factory<Matrix, TrivialGroup> {
        static typename model_traits<Matrix, TrivialGroup>::model_ptr parse
        (Lattice const & lattice, BaseParameters & model)
        {
            if (model.get<std::string>("MODEL") == std::string("optical_lattice"))
                return typename model_traits<Matrix, TrivialGroup>::model_ptr(
                            new OpticalLatticeNull<Matrix>(lattice, model)
                       );
            else {
                throw std::runtime_error("Don't know this model!");
                return typename model_traits<Matrix, TrivialGroup>::model_ptr();
            }
        }
    };
    
}