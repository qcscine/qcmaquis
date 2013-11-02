/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/chem/model_qc.h"

template<class Matrix>
struct model_factory<Matrix, TwoU1PG> {
    static typename model_traits<Matrix, TwoU1PG>::model_ptr parse
    (Lattice const & lattice, BaseParameters & model)
    {
        if (model["MODEL"] == std::string("quantum_chemistry"))
            return typename model_traits<Matrix, TwoU1PG>::model_ptr(
                    new qc_model<Matrix, TwoU1PG>(lattice, model)
                   );

        else {
            throw std::runtime_error("Don't know this model!");
            return typename model_traits<Matrix, TwoU1PG>::model_ptr();
        }
    }
};
