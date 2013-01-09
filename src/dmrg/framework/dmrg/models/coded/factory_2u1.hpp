/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/coded/models_2u1.hpp"
#include "dmrg/models/chem/model_qc.h"

template<class Matrix>
struct model_factory<Matrix, TwoU1> {
    static typename model_traits<Matrix, TwoU1>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        if (model.get<std::string>("MODEL") == std::string("fermion Hubbard"))
            return typename model_traits<Matrix, TwoU1>::model_ptr(
                        new FermiHubbardTwoU1<Matrix>(lattice, model)
                   );
        else if (model.get<std::string>("MODEL") == std::string("quantum_chemistry"))
            return typename model_traits<Matrix, TwoU1>::model_ptr(
                    new qc_model<Matrix>(lattice, model)
                   );

        else {
            throw std::runtime_error("Don't know this model!");
            return typename model_traits<Matrix, TwoU1>::model_ptr();
        }
    }
};
