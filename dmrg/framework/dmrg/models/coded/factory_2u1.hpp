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
struct coded_model_factory<Matrix, TwoU1> {
    static boost::shared_ptr<model_impl<Matrix, TwoU1> > parse
    (Lattice const & lattice, BaseParameters & model)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TwoU1> > impl_ptr;
        if (model["MODEL"] == std::string("fermion Hubbard"))
            return impl_ptr( new FermiHubbardTwoU1<Matrix>(lattice, model) );
        else
            if (model["MODEL"] == std::string("quantum_chemistry"))
            return impl_ptr( new qc_model<Matrix, TwoU1>(lattice, model) );

        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
