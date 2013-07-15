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
struct cont_model_factory<Matrix, TwoU1> {
    static typename model_traits<Matrix, TwoU1>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        std::string model_str = model.is_set("model") ? "model" : "MODEL";
        if (model[model_str] == std::string("fermi_optical_lattice"))
            return typename model_traits<Matrix, TwoU1>::model_ptr(
                        new FermiOpticalLattice<Matrix>(lattice, model)
                   );
        else if (model[model_str] == std::string("optical_lattice_cons_dm"))
            return typename model_traits<Matrix, TwoU1>::model_ptr( new DMOpticalLatticeTwoU1<Matrix>(lattice, model) );
        else {
            throw std::runtime_error("Don't know this model!");
            return typename model_traits<Matrix, TwoU1>::model_ptr();
        }
    }
};
