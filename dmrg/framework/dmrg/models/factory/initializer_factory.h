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
//    template <class Matrix, class SymmGroup>
//    struct call_linear_init {
//        static typename Model<Matrix,SymmGroup>::initializer_ptr call()
//        {
//            throw std::runtime_error("Linear MPS init is available only for U1 symmetry group.");
//            BaseParameters bp;
//            return typename Model<Matrix,SymmGroup>::initializer_ptr(new default_mps_init<Matrix, SymmGroup>(bp));
//        }
//    };
//    template <class Matrix>
//    struct call_linear_init<Matrix, U1> {
//        static typename Model<Matrix,U1>::initializer_ptr call()
//        {
//            return typename Model<Matrix,U1>::initializer_ptr(new linear_mps_init<Matrix>());
//        }
//    };
    
    template <class Matrix, class SymmGroup>
    struct call_hf_init {
        static typename Model<Matrix,SymmGroup>::initializer_ptr call(BaseParameters parms, BaseParameters model,
                                                                      std::vector<Index<SymmGroup> > const& phys_dims,
                                                                      typename SymmGroup::charge right_end,
                                                                      std::vector<int> const& site_type)
        {
            throw std::runtime_error("HF MPS init is available only for TwoU1 or TwoU1PG symmetry group.");
            return typename Model<Matrix,SymmGroup>::initializer_ptr(new default_mps_init<Matrix, SymmGroup>(parms, model, phys_dims, right_end, site_type));
        }
    };
    template <class Matrix>
    struct call_hf_init<Matrix, TwoU1> {
        static typename Model<Matrix,TwoU1>::initializer_ptr call(BaseParameters parms, BaseParameters model,
                                                                  std::vector<Index<TwoU1> > const& phys_dims,
                                                                  TwoU1::charge right_end,
                                                                  std::vector<int> const& site_type)
        {
            return typename Model<Matrix,TwoU1>::initializer_ptr(new hf_mps_init<Matrix, TwoU1>(parms, model, phys_dims, right_end, site_type));
        }
    };
    template <class Matrix>
    struct call_hf_init<Matrix, TwoU1PG> {
        static typename Model<Matrix,TwoU1PG>::initializer_ptr call(BaseParameters parms, BaseParameters model,
                                                                    std::vector<Index<TwoU1PG> > const& phys_dims,
                                                                    TwoU1PG::charge right_end,
                                                                    std::vector<int> const& site_type)
        {
            return typename Model<Matrix,TwoU1PG>::initializer_ptr(new hf_mps_init<Matrix, TwoU1PG>(parms, model, phys_dims, right_end, site_type));
        }
    };
}

template <class Matrix, class SymmGroup>
typename model_impl<Matrix,SymmGroup>::initializer_ptr
model_impl<Matrix,SymmGroup>::initializer(Lattice const& lat, BaseParameters & parms, BaseParameters & model) const
{
    typename SymmGroup::charge initc = this->total_quantum_numbers(model);
    
    int max_site_type = 0;
    std::vector<int> site_types(lat.size(), 0);
    for (int p = 0; p < lat.size(); ++p) {
        site_types[p] = lat.get_prop<int>("type", p);
        max_site_type = std::max(site_types[p], max_site_type);
    }
    
    std::cout << "site_types: ";
    std::copy(site_types.begin(), site_types.end(), std::ostream_iterator<int>(std::cout, " "));
    std::cout << std::endl;
    
    std::vector<Index<SymmGroup> > site_bases(max_site_type+1);
    for (int type = 0; type < site_bases.size(); ++type) {
        site_bases[type] = this->phys_dim(type);
        std::cout << "phys["<<type <<"]: " << site_bases[type] << std::endl;
    }
    
    if (parms["init_state"] == "default")
        return initializer_ptr(new default_mps_init<Matrix, SymmGroup>(parms, model, site_bases, initc, site_types));
    
//    else if (params["init_state"] == "linear")
//        return detail::call_linear_init<Matrix, SymmGroup>::call();
    
    else if (parms["init_state"] == "const")
        return initializer_ptr(new const_mps_init<Matrix, SymmGroup>(parms, model, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "thin")
        return initializer_ptr(new thin_mps_init<Matrix, SymmGroup>(parms, model, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "thin_const")
        return initializer_ptr(new thin_const_mps_init<Matrix, SymmGroup>(parms, model, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "basis_state")
        return initializer_ptr(new basis_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "basis_state_generic")
        return initializer_ptr(new basis_mps_init_generic<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    
    else if (parms["init_state"] == "coherent")
        return initializer_ptr(new coherent_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "basis_state_dm")
        return initializer_ptr(new basis_dm_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "coherent_dm")
        return initializer_ptr(new coherent_dm_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    
    else if (parms["init_state"] == "hf")
        return detail::call_hf_init<Matrix, SymmGroup>::call(parms, model, site_bases, initc, site_types);
    
    else {
        throw std::runtime_error("Don't know this initial state.");
        return initializer_ptr();
    }

}

#endif
