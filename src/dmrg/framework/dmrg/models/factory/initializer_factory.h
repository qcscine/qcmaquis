/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H
#define MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H

#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/models/chem/mps_init_hf.hpp"

namespace detail {
    template <class Matrix, class SymmGroup>
    struct call_linear_init {
        static typename Model<Matrix,SymmGroup>::initializer_ptr call()
        {
            throw std::runtime_error("Linear MPS init is available only for U1 symmetry group.");
            return typename Model<Matrix,SymmGroup>::initializer_ptr(new default_mps_init<Matrix, SymmGroup>());
        }
    };
    template <class Matrix>
    struct call_linear_init<Matrix, U1> {
        static typename Model<Matrix,U1>::initializer_ptr call()
        {
            return typename Model<Matrix,U1>::initializer_ptr(new linear_mps_init<Matrix>());
        }
    };
    
    template <class Matrix, class SymmGroup>
    struct call_hf_init {
        static typename Model<Matrix,SymmGroup>::initializer_ptr call(BaseParameters & params)
        {
            throw std::runtime_error("HF MPS init is available only for TwoU1 symmetry group.");
            return typename Model<Matrix,SymmGroup>::initializer_ptr(new default_mps_init<Matrix, SymmGroup>());
        }
    };
    template <class Matrix>
    struct call_hf_init<Matrix, TwoU1> {
        static typename Model<Matrix,TwoU1>::initializer_ptr call(BaseParameters & params)
        {
            return typename Model<Matrix,TwoU1>::initializer_ptr(new hf_mps_init<Matrix>(params));
        }
    };
}

template <class Matrix, class SymmGroup>
typename Model<Matrix,SymmGroup>::initializer_ptr Model<Matrix,SymmGroup>::initializer(BaseParameters & params) const
{
    if (params.get<std::string>("init_state") == "default")
        return initializer_ptr(new default_mps_init<Matrix, SymmGroup>());
    
    else if (params.get<std::string>("init_state") == "linear")
        return detail::call_linear_init<Matrix, SymmGroup>::call();
    
    else if (params.get<std::string>("init_state") == "const")
        return initializer_ptr(new const_mps_init<Matrix, SymmGroup>());
    
    else if (params.get<std::string>("init_state") == "thin")
        return initializer_ptr(new thin_mps_init<Matrix, SymmGroup>());
    
    else if (params.get<std::string>("init_state") == "thin_const")
        return initializer_ptr(new thin_const_mps_init<Matrix, SymmGroup>());
    
    else if (params.get<std::string>("init_state") == "basis_state")
        return initializer_ptr(new basis_mps_init<Matrix, SymmGroup>(params));
    
    else if (params.get<std::string>("init_state") == "basis_state_generic")
        return initializer_ptr(new basis_mps_init_generic<Matrix, SymmGroup>(params));
    
    else if (params.get<std::string>("init_state") == "coherent")
        return initializer_ptr(new coherent_mps_init<Matrix, SymmGroup>(params));
    
    else if (params.get<std::string>("init_state") == "basis_state_dm")
        return initializer_ptr(new basis_dm_mps_init<Matrix, SymmGroup>(params));
    
    else if (params.get<std::string>("init_state") == "coherent_dm")
        return initializer_ptr(new coherent_dm_mps_init<Matrix, SymmGroup>(params));
    
    else if (params.get<std::string>("init_state") == "hf")
        return detail::call_hf_init<Matrix, SymmGroup>::call(params);
    
    else {
        throw std::runtime_error("Don't know this initial state.");
        return initializer_ptr();
    }

}

#endif
