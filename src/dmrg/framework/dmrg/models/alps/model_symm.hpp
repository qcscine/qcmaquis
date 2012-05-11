/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_ALPS_MODEL_SYMM_H
#define APP_ALPS_MODEL_SYMM_H

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

#include "dmrg/block_matrix/symmetry.h"


namespace detail {
	template <class T>
	int to_integer (alps::half_integer<T> const & qn_value)
	{
		return (qn_value.get_twice() % 2 == 0) ? alps::to_integer(qn_value) : qn_value.get_twice();
	}
}

template <class SymmGroup>
typename SymmGroup::charge convert_alps (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn);

template <class SymmGroup>
typename SymmGroup::charge init_charge (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn);

template <class SymmGroup>
std::map<typename SymmGroup::charge,std::size_t> init_qn_charges
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states);

template <class SymmGroup>
std::map<alps::site_state<short>, std::pair<typename SymmGroup::charge, std::size_t> > init_coords
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states);

#endif
