/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#ifndef MPS_INIT_HF_HPP
#define MPS_INIT_HF_HPP

#include "dmrg/mp_tensors/compression.h"
#include <boost/type_traits.hpp>

template<class Matrix, class SymmGroup>
class hf_mps_init : public mps_initializer<Matrix, SymmGroup>
{
public:
    typedef Lattice::pos_t pos_t;
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;
    typedef typename mps_initializer<Matrix, SymmGroup>::MPSVector MPSVector;
    typedef std::set<charge> container_type;

    hf_mps_init(BaseParameters parms_,
                std::vector<Index<SymmGroup> > const& phys_dims_,
                charge right_end,
                std::vector<int> const& site_type)
    : parms(parms_)
    , phys_dims(phys_dims_)
    , site_types(site_type)
    , di(parms, phys_dims_, right_end, site_type)
    {}

    BaseParameters parms;
    std::vector<Index<SymmGroup> > phys_dims;
    std::vector<int> site_types;
    default_mps_init<Matrix, SymmGroup> di;

    void operator()(MPS<Matrix, SymmGroup> & mps)
    {
        di.init_sectors(mps, 5, false, 1.);

        std::vector<std::size_t> hf_init = parms["hf_occ"];

        init_hf_mps(mps, hf_init, symm_traits::HasSU2<SymmGroup>());
    }

    void operator()(MPSVector & mps_vector)
    {
      // HF guesses for the SA case should be separated with the semicolons
      std::string hf_init_all = parms["hf_occ"];
      std::vector<std::string> hf_init_sa;

      boost::split(hf_init_sa, hf_init_all, boost::is_any_of(";"));

      // number of HF guesses must be the same as the number of the SA states

      if (hf_init_sa.size() != mps_vector.size())
        throw std::runtime_error("Number of HF guesses does not match the number of states!");

      for (typename MPSVector::iterator it = mps_vector.begin(); it != mps_vector.end(); ++it)
      {
        // parameters::proxy takes care of the conversion from string of integers to std::vector<std::size_t>
        parameters::proxy parm(*(hf_init_sa.begin() + std::distance(mps_vector.begin(), it)));
        std::vector<std::size_t> hf_init = parm.as<std::vector<std::size_t> >();

        di.init_sectors(*it, 5, false, 1.);

        maquis::cout << "Initializing state " << std::distance(mps_vector.begin(), it) << " with HF occupation vector ";
        std::copy(hf_init.begin(), hf_init.end(), std::ostream_iterator<std::size_t>(maquis::cout, ""));
        maquis::cout << std::endl;

        init_hf_mps(*it, hf_init, symm_traits::HasSU2<SymmGroup>());

      }
    }

    // non-SU2 version of HF initialiser
    void init_hf_mps(MPS<Matrix,SymmGroup> & mps, std::vector<std::size_t> hf_init, boost::false_type)
    {
      std::vector<pos_t> order(mps.length());
      if (!parms.is_set("orbital_order"))
          for (pos_t p = 0; p < mps.length(); ++p)
              order[p] = p+1;
      else
          order = parms["orbital_order"].template as<std::vector<pos_t> >();

      std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);

      if (hf_init.size() != mps.length())
          throw std::runtime_error("HF occupation vector length != MPS length\n");

      charge max_charge = SymmGroup::IdentityCharge;
      for (pos_t i = 0; i < mps.length(); ++i)
      {
          mps[i].multiply_by_scalar(0.0);

          size_t sc_input = hf_init[order[i]];
          typename SymmGroup::charge site_charge(0);

          size_t loc_dim = phys_dims[0].size();

          if (sc_input > loc_dim)
              throw std::runtime_error(
                  "The hf occ exceeds local basis dimension\n"
              );

          site_charge = phys_dims[site_types[i]][loc_dim-sc_input].first;

          max_charge = SymmGroup::fuse(max_charge, site_charge);

          #ifndef NDEBUG
          maquis::cout << "site " << i << " activating sector " << max_charge << std::endl;
          #endif

          // Set largest charge sector = all 1
          size_t max_pos = mps[i].data().left_basis().position(max_charge);
          if (max_pos >= mps[i].data().n_blocks()) {
              maquis::cout << "ERROR: Symmetry block " << max_charge << " not found\n";
              maquis::cout << "site " << i << ", site_charge " << site_charge << ", cumulated_charge "
                        << max_charge << ", block_pos: "
                        << max_pos << ", number of blocks: " << mps[i].data().n_blocks() << std::endl;

              maquis::cout << "This error occurs if the specified HF determinant is not in the same symmetry sector as the target state\n";
              exit(1);
          }
          Matrix & mfirst = mps[i].data()[max_pos];
          size_t nrow = mfirst.num_rows();
          size_t ncol = mfirst.num_cols();
          mps[i].data()[max_pos] = Matrix(nrow, ncol, 1.);

          mps[i].multiply_by_scalar(1. / mps[i].scalar_norm());
      }

      //mps = compression::l2r_compress(mps, init_bond_dimension, 1e-6);

      //maquis::cout << "\nMPS AFTER COMPRESSION:\n";
      //for(int i = 0; i < mps.length(); ++i) {
      //    maquis::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
      //    maquis::cout << mps[i].scalar_norm() << std::endl;
      //}
    }

    // SU2 version of HF initialiser
    void init_hf_mps(MPS<Matrix,SymmGroup> & mps, std::vector<std::size_t> hf_init, boost::true_type)
    {
        std::vector<pos_t> order(mps.length());
        if (!parms.is_set("orbital_order"))
            for (pos_t p = 0; p < mps.length(); ++p)
                order[p] = p+1;
        else
            order = parms["orbital_order"].template as<std::vector<pos_t> >();

        std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);

        if (hf_init.size() != mps.length())
            throw std::runtime_error("HF occupation vector length != MPS length\n");

        container_type bond_charges;
        bond_charges.insert(SymmGroup::IdentityCharge);
        for (pos_t i = 0; i < mps.length(); ++i)
        {
            mps[i].multiply_by_scalar(0.0);

            size_t sc_input = hf_init[order[i]];
            container_type site_charges;

            if (sc_input > 4)
                throw std::runtime_error(
                    "The hf_occ format has been changed to: 1=empty, 2=down, 3=up, 4=updown\n (not cumulative anymore)\n"
                );
            if (phys_dims[site_types[i]].size() != 4) throw std::runtime_error("HF init expects 4 states per orbital\n");

            switch(sc_input) {
                case 4:
                    site_charges.insert(phys_dims[site_types[i]][0].first); // doubly-occ
                    break;
                case 3:
                    site_charges.insert(phys_dims[site_types[i]][1].first); // singly-occ
                    site_charges.insert(phys_dims[site_types[i]][2].first); // singly-occ
                    break;
                case 2:
                    site_charges.insert(phys_dims[site_types[i]][1].first); // singly-occ
                    site_charges.insert(phys_dims[site_types[i]][2].first); // singly-occ
                    break;
                case 1:
                    site_charges.insert(phys_dims[site_types[i]][3].first); // empty
                    break;
            }

            container_type next_bond_charges;
            for (typename container_type::const_iterator it = bond_charges.begin(); it != bond_charges.end(); ++it)
                for (typename container_type::const_iterator it2 = site_charges.begin(); it2 != site_charges.end(); ++it2)
                {
                    charge sector = SymmGroup::fuse(*it, *it2);
                    if (mps[i].col_dim().has(sector))
                        next_bond_charges.insert(sector);
                }

            #ifndef NDEBUG
            for (typename container_type::const_iterator it = next_bond_charges.begin(); it != next_bond_charges.end(); ++it)
                maquis::cout << "site " << i << " activating sector " << *it << std::endl;

            maquis::cout << std::endl;
            #endif

            for (typename container_type::const_iterator it = next_bond_charges.begin(); it != next_bond_charges.end(); ++it)
            {
                charge max_charge = *it;
                size_t max_pos = mps[i].data().left_basis().position(max_charge);
                if (max_pos >= mps[i].data().n_blocks()) {
                    maquis::cout << "ERROR: Symmetry block " << max_charge << " not found\n";
                    maquis::cout << "site " << i << ", cumulated_charge "
                             << max_charge << ", block_pos: "
                             << max_pos << ", number of blocks: " << mps[i].data().n_blocks() << std::endl;

                    maquis::cout << "This error occurs if the specified HF determinant is not in the same symmetry sector as the target state\n";
                    exit(1);
                }
                Matrix & mfirst = mps[i].data()[max_pos];
                size_t nrow = mfirst.num_rows();
                size_t ncol = mfirst.num_cols();
                // Set current charge sector = all 1
                mps[i].data()[max_pos] = Matrix(nrow, ncol, 1.);
            }

            mps[i].multiply_by_scalar(1. / mps[i].scalar_norm());

#ifndef NDEBUG
            maquis::cout << "mps[" << i << "]:\n" << mps[i] << std::endl;
#endif
            std::swap(next_bond_charges, bond_charges);
        }
    }
};

#endif
