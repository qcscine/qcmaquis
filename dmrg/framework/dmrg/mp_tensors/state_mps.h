/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_STATE_MPS_H
#define MAQUIS_DMRG_STATE_MPS_H

#include "dmrg/mp_tensors/mps.h"
#include "mps_sectors.h"
#include <boost/tuple/tuple.hpp>

template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> state_mps(std::vector<boost::tuple<typename SymmGroup::charge, int> > const & state,
                                 std::vector<Index<SymmGroup> > const& phys_dims,
                                 std::vector<int> const& site_type,
                                 typename SymmGroup::charge right_end = SymmGroup::IdentityCharge,
                                 int mdim=1)
{
    // Types and variable definition
    typedef typename SymmGroup::charge charge;
    MPS<Matrix, SymmGroup> mps(state.size());
    Index<SymmGroup> curr_i;
    std::vector< Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, mdim);
    // -- MAIN LOOP --
    // The overall structure of the algorithm is the following:
    // One first generates the index (IdentityCharge,mdim), where mdim is the number of renormalized
    // block states (1 in the default case). This is the index with which we start from the left.
    // For a given left index and a given physical basis state, the "acceptable" QN for the right
    // renormalized basis are univocally determined by the product of the two set of QNs.
    // Therefore, since we have a single ONV as a starting point, we will have
    // all 1x1 tensors, but the physical dimensions will be in general != from 0.
    curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, mdim));
    for (int i = 0; i < state.size(); ++i)
    {
        // Computes the symmetry block of the next dimension and HARDCODED its value to 1
        charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
        Index<SymmGroup> new_i;
        new_i.insert(std::make_pair(newc, mdim));
        // Get the product basis between the physical basis and the symmetry block of the left renormalized basis
        ProductBasis<SymmGroup> left(phys_dims[site_type[i]], allowed[i]);
        mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], false, 0);
        // Finds out where to put the 1.0 in the MPS. Retrieve, from the ProductBasis object, how the row index was
        // decomposed in terms of left auxiliary basis and physical basis.
        size_t b_in = left(boost::get<0>(state[i]), curr_i[0].first) + boost::get<1>(state[i]) * curr_i[0].second;
        assert (allowed[i+1].has(newc));
        size_t b_out = 0;
        mps[i].make_left_paired();
        // Populates the MPS
        block_matrix<Matrix, SymmGroup> & block = mps[i].data();
        Matrix &m = block(newc, new_i[0].first);
        m(b_in, b_out) = 1.;
        curr_i = new_i;
    }
    return mps;
}

/** @brief Same as above, but does not populate only the i-th position in the ONV, but all positions up to i */
template <class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> state_mps_const(std::vector<boost::tuple<typename SymmGroup::charge, int> > const & state,
                                       std::vector<Index<SymmGroup> > const& phys_dims, std::vector<int> const& site_type,
                                       typename SymmGroup::charge right_end = SymmGroup::IdentityCharge, bool fillRand=false,
                                       int mMax=1)
{
  int mdim=1;
  // Types definition
  typedef typename SymmGroup::charge charge;
  MPS<Matrix, SymmGroup> mps(state.size());
  Index<SymmGroup> curr_i;
  std::vector< Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, mdim);
  curr_i.insert(std::make_pair(SymmGroup::IdentityCharge, mdim));
  for (int i = 0; i < state.size(); ++i) {
    // Computes the symmetry block of the next dimension and HARDCODED its value to 1
    charge newc = SymmGroup::fuse(curr_i[0].first, boost::get<0>(state[i]));
    assert (allowed[i+1].has(newc));
    Index<SymmGroup> new_i;
    new_i.insert(std::make_pair(newc, mdim));
    // Get the product basis between the physical basis and the symmetry block of the left renormalized basis
    ProductBasis<SymmGroup> left(phys_dims[site_type[i]], allowed[i]);
    mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], false, 0.);
    mps[i].make_left_paired();
    // Finds out where to put the 1.0 in the MPS.
    auto refIndex = boost::get<1>(state[i]);
    int lowerBound, upperBound;
    if (refIndex < 0) {
      lowerBound = -refIndex;
      upperBound = -refIndex;
    }
    else {
      lowerBound = 0;
      upperBound = refIndex;
    }
    for (int iBlock = lowerBound; iBlock <= upperBound; iBlock++) {
      auto b_in = left(boost::get<0>(state[i]), curr_i[0].first) + iBlock * curr_i[0].second;
      mps[i].data()(newc, new_i[0].first)(b_in, 0) = fillRand ? dmrg_random::uniform(0., 1.) : 1.;
    }
    curr_i = new_i;
  }
  return mps;
}

// Special case NU1
// @brief Same as above, but does not populate only the i-th position in the ONV, but all positions up to i for each mode
template <class Matrix, int N>
MPS<Matrix, NU1_template<N>> state_mps_const(std::vector<boost::tuple<typename NU1_template<N>::charge, int> > const & state,
                                       std::vector<Index<NU1_template<N>> > const& phys_dims, std::vector<int> const& site_type,
                                       typename NU1_template<N>::charge right_end = NU1_template<N>::IdentityCharge, bool fillRand=false,
                                       int mMax = 1)
{
  // Types definition
  using SymmGroup = NU1_template<N>;
  typedef typename SymmGroup::charge charge;
  MPS<Matrix, SymmGroup> mps(state.size());
  std::vector< Index<SymmGroup> > allowed = allowed_sectors(site_type, phys_dims, right_end, mMax);
  // loop over all sites
  for (int i = 0; i < state.size(); ++i) {
    auto populate = boost::get<1>(state[i]); // state should be now pair of charge and int as bool (const/rand or zero)
    if (populate) {
      mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], fillRand, 1.0);
    } else {
      mps[i] = MPSTensor<Matrix, SymmGroup>(phys_dims[site_type[i]], allowed[i], allowed[i+1], false, 0);
      mps[i].make_left_paired();
      // Get the product basis between the physical basis and the symmetry block of the left renormalized basis
      ProductBasis<SymmGroup> left(phys_dims[site_type[i]], allowed[i]);
      for (int iBlock = 0; iBlock < allowed[i].size(); iBlock++) {
        charge newc = allowed[i][iBlock].first;
        if (allowed[i+1].has(newc)){
          auto offset = left(SymmGroup::IdentityCharge, newc);
          block_matrix<Matrix, SymmGroup> & block = mps[i].data();
          Matrix &m = block(newc, newc);
          for (int rowInBlock = offset; rowInBlock < allowed[i][iBlock].second + offset; rowInBlock++){
            for (int columnInBlock = 0; columnInBlock < allowed[i+1].size_of_block(newc); columnInBlock++) {
              m(rowInBlock, columnInBlock) = fillRand ? dmrg_random::uniform(0., 1.) : 1.;
            }
          }
        }
      }
    }
  }
  return mps;
}

#endif
