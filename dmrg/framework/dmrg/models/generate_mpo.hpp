/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "dmrg/models/generate_mpo/mpo_maker.hpp"
#include "dmrg/models/generate_mpo/tagged_mpo_maker_optim.hpp"
#include "dmrg/models/generate_mpo/corr_maker.hpp"
#include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"


template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_mpo(Lattice const& lat, Model<Matrix, SymmGroup> & model)
{
    model.create_terms();
    generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpom(lat, model);
    MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
    return mpo;
}

#endif
