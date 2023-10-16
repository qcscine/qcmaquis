/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H
#define MAQUIS_DMRG_MODELS_INITIALIZER_FACTORY_H

#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/models/chem/mps_init_hf.hpp"

namespace detail {

template <class Matrix, class SymmGroup, class = void>
struct call_hf_init {
    static typename Model<Matrix,SymmGroup>::initializer_ptr call(BaseParameters parms,
                                                                  std::vector<Index<SymmGroup> > const& phys_dims,
                                                                  typename SymmGroup::charge right_end,
                                                                  std::vector<int> const& site_type)
    {
        throw std::runtime_error("No HF MPS init available for this symmetry group.");
        return typename Model<Matrix,SymmGroup>::initializer_ptr(new default_mps_init<Matrix, SymmGroup>(parms, phys_dims, right_end, site_type));
    }
};

template <class Matrix, class SymmGroup>
struct call_hf_init<Matrix, SymmGroup, symm_traits::enable_if_chemmodel_t<SymmGroup> > {
    static typename Model<Matrix,SymmGroup>::initializer_ptr call(BaseParameters parms,
                                                                  std::vector<Index<SymmGroup> > const& phys_dims,
                                                                  typename SymmGroup::charge right_end,
                                                                  std::vector<int> const& site_type)
    {
        return typename Model<Matrix,SymmGroup>::initializer_ptr(new hf_mps_init<Matrix, SymmGroup>(parms, phys_dims, right_end, site_type));
    }
};

} // namespace detail

template <class Matrix, class SymmGroup>
typename model_impl<Matrix,SymmGroup>::initializer_ptr
model_impl<Matrix,SymmGroup>::initializer(Lattice const& lat, BaseParameters & parms) const
{
    typename SymmGroup::charge initc = this->total_quantum_numbers(parms);

    int max_site_type = 0;
    std::vector<int> site_types(lat.size(), 0);
    for (int p = 0; p < lat.size(); ++p) {
        site_types[p] = lat.get_prop<int>("type", p);
        max_site_type = std::max(site_types[p], max_site_type);
    }

    // maquis::cout << "site_types: ";
    // std::copy(site_types.begin(), site_types.end(), maquis::ostream_iterator<int>(maquis::cout, " "));
    // maquis::cout << std::endl;

    std::vector<Index<SymmGroup> > site_bases(max_site_type+1);
    for (int type = 0; type < site_bases.size(); ++type) {
        site_bases[type] = this->phys_dim(type);
        // maquis::cout << "phys["<< type <<"]: " << site_bases[type] << std::endl;
    }

    // Generation of the initializer
    if (parms["init_type"] == "default")
        return initializer_ptr(new default_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    else if (parms["init_type"] == "const")
        return initializer_ptr(new const_mps_init<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    else if (parms["init_type"] == "basis_state")
        return initializer_ptr(new basis_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    else if (parms["init_type"] == "basis_state_generic")
        return initializer_ptr(new basis_mps_init_generic<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    else if (parms["init_type"] == "basis_state_generic_const")
        return initializer_ptr(new basis_mps_init_generic_const<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    else if (parms["init_type"] == "basis_state_generic_default")
        return initializer_ptr(new basis_mps_init_generic_default<Matrix, SymmGroup>(parms, site_bases, initc, site_types));
    else if (parms["init_type"] == "coherent")
        return initializer_ptr(new coherent_mps_init<Matrix, SymmGroup>(parms, site_bases, site_types));
    else if (parms["init_type"] == "hf")
        return detail::call_hf_init<Matrix, SymmGroup>::call(parms, site_bases, initc, site_types);
    else {
        throw std::runtime_error("Don't know this initial state.");
        return initializer_ptr();
    }
}

#endif
