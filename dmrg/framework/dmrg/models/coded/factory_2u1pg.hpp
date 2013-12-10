/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/models/chem/model_qc.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1PG> {
    static boost::shared_ptr<model_impl<Matrix, TwoU1PG> > parse
    (Lattice const & lattice, BaseParameters & model)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TwoU1PG> > impl_ptr;
        if (model["MODEL"] == std::string("quantum_chemistry")) {
            if (model["LATTICE"] == std::string("quantum_chemistry"))
                throw std::runtime_error("Please use \"LATTICE = orbitals\" for quantum_chemistry\n");

            return impl_ptr( new qc_model<Matrix, TwoU1PG>(lattice, model) );
        }

        else {
            throw std::runtime_error("Don't know this model!\n");
            return impl_ptr();
        }
    }
};
