/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_MODELS_H
#define MAQUIS_DMRG_MODELS_MODELS_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/hamiltonian.h"
#include "dmrg/models/measurements.h"

#include <boost/shared_ptr.hpp>


namespace app {
    
    template <class SymmGroup>
    typename SymmGroup::charge init_qn (BaseParameters &);
    
    template <class Matrix, class SymmGroup>
    class Model {
    public:
        virtual Hamiltonian<Matrix, SymmGroup> H () const=0;
        virtual Measurements<Matrix, SymmGroup> measurements () const=0;
        virtual typename SymmGroup::charge initc(BaseParameters & parms)
        {
            return init_qn<SymmGroup>(parms);
        }
    };
    
    
    
    template <class Matrix, class SymmGroup>
    struct model_traits {
        typedef Model<Matrix, SymmGroup> model;
        typedef boost::shared_ptr<model> model_ptr;
        
    };
    
}

#endif
