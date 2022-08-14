/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
 *               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
 * Copyright (C) 2018 Institute for Theoretical Chemistry, ETH Zurich
 *               2018-2019    Stefan Knecht <stknecht@ethz.ch>
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

#ifndef MEASUREMENTS_TAGGED_NRANKRDM_REL_H
#define MEASUREMENTS_TAGGED_NRANKRDM_REL_H

#include <algorithm>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>
#include "dmrg/block_matrix/symmetry/nu1pg.h"
#include "dmrg/models/measurement.h"
#include "dmrg/utils/checks.h"
#include "dmrg/models/MolecularHamiltonians/su2u1/term_maker.h"
#include "dmrg/models/MolecularHamiltonians/transform_symmetry.hpp"
#include "measurements_details.h"

namespace measurements {

template <class Matrix, class SymmGroup>
class TaggedNRankRDM<Matrix, SymmGroup, symm_traits::enable_if_u1dg_t<SymmGroup> >
    : public measurement<Matrix, SymmGroup>
{
  typedef measurement<Matrix, SymmGroup> base;
  typedef typename Model<Matrix, SymmGroup>::term_descriptor term_descriptor;
  typedef Lattice::pos_t pos_t;
  typedef std::vector<pos_t> positions_type;
  typedef typename base::op_t op_t;
  typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
  typedef typename Matrix::value_type value_type;
  typedef std::vector<tag_type> tag_vec;
  typedef std::vector<tag_vec> bond_term;
  typedef std::pair<std::vector<tag_vec>, value_type> scaled_bond_term;

public:
  TaggedNRankRDM(std::string const& name_, const Lattice & lat, std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                 tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                 bool half_only_, positions_type const& positions_ = positions_type(), std::string const& ckp_ = std::string(""))
    : base(name_), lattice(lat), tag_handler(tag_handler_), positions_first(positions_), identities(identities_),
      fillings(fillings_), operator_terms(ops_), bra_ckp(ckp_)
  {
    pos_t extent = lattice.size();
    if (positions_first.size() == 0)
    {
        positions_first.resize(extent);
        std::iota(positions_first.begin(), positions_first.end(), 0);
    }
    this->cast_to_real = false;
  }

  void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
  {
    this->vector_results.clear();
    this->labels.clear();
    this->labels_num.clear();
    MPS<Matrix, SymmGroup> bra_mps;
    if (bra_ckp != "") {
        if(boost::filesystem::exists(bra_ckp))
            load(bra_ckp, bra_mps);
        else
            throw std::runtime_error("The bra checkpoint file " + bra_ckp + " was not found\n");
    }

    maquis::cout << " measuring in rel version of tagged_nrank " << std::endl;

    if (operator_terms[0].first.size() == 2)
      measure_correlation(bra_mps, ket_mps);
    else if (operator_terms[0].first.size() == 4)
      measure_2rdm(bra_mps, ket_mps);
    //else if (operator_terms[0].first.size() == 6)
    //    measure_3rdm(bra_mps, ket_mps);
    //else if (operator_terms[0].first.size() == 8)
    //    measure_4rdm(bra_mps, ket_mps);
    else
        throw std::runtime_error("relativistic correlation measurements at the moment supported with 2 and 4 operators, size is "
                                  + boost::lexical_cast<std::string>(operator_terms[0].first.size()));
  }

protected:

  measurement<Matrix, SymmGroup>* do_clone() const { return new TaggedNRankRDM(*this); }

  void measure_correlation(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                           MPS<Matrix, SymmGroup> const & ket_mps)
  {
    // Test if a separate bra state has been specified
    bool bra_neq_ket = (dummy_bra_mps.length() > 0);
    //MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
    MPS<Matrix, SymmGroup> bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
    MPS<Matrix, SymmGroup> ket_mps_local = ket_mps;

    #ifdef MAQUIS_OPENMP
    #pragma omp parallel for schedule(dynamic) firstprivate(ket_mps_local, bra_mps)
    #endif
    for (std::size_t i = 0; i < positions_first.size(); ++i) {
      pos_t p1 = positions_first[i];
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));

      std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
      std::vector<std::vector<pos_t> > num_labels;
      for (pos_t p2 = p1; p2 < lattice.size(); ++p2)
      {
        std::vector<pos_t> positions = {p1, p2};

        tag_vec operators(2);
        operators[0] = operator_terms[0].first[0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
        operators[1] = operator_terms[0].first[1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];

        // check if term is allowed by symmetry
        term_descriptor term = generate_mpo::arrange_operators(positions, operators, tag_handler_local);

        //if(measurements_details::checkpg<SymmGroup>()(term, tag_handler_local, lattice))
        {
          MPO<Matrix, SymmGroup> mpo = generate_mpo::sign_and_fill(term, identities, fillings, tag_handler_local, lattice);
          typename MPS<Matrix, SymmGroup>::scalar_type value = operator_terms[0].second * expval(bra_mps, ket_mps_local, mpo);
          dct.push_back(value);
          num_labels.push_back(order_labels(lattice, positions));
        }
        //else {
        //    dct.push_back(0.0);
        //    num_labels.push_back(positions);
        //}
      }

      std::vector<std::string> lbt = label_strings(num_labels);

      // save results and labels
      #ifdef MAQUIS_OPENMP
      #pragma omp critical
      #endif
      {
        this->vector_results.reserve(this->vector_results.size() + dct.size());
        std::copy(dct.rbegin(), dct.rend(), std::back_inserter(this->vector_results));

        this->labels.reserve(this->labels.size() + dct.size());
        std::copy(lbt.rbegin(), lbt.rend(), std::back_inserter(this->labels));

        this->labels_num.reserve(this->labels_num.size() + dct.size());
        std::copy(num_labels.rbegin(), num_labels.rend(), std::back_inserter(this->labels_num));
      }
    }
  }

  void measure_2rdm(MPS<Matrix, SymmGroup> const & dummy_bra_mps, MPS<Matrix, SymmGroup> const & ket_mps)
  {
    // Test if a separate bra state has been specified
    bool bra_neq_ket = (dummy_bra_mps.length() > 0);
    //MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
    MPS<Matrix, SymmGroup> bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
    MPS<Matrix, SymmGroup> ket_mps_local = ket_mps;
    #ifdef MAQUIS_OPENMP
    #pragma omp parallel for collapse(1) schedule(dynamic) firstprivate(ket_mps_local, bra_mps)
    #endif
    for (pos_t p1 = 0; p1 < lattice.size(); ++p1)
    for (pos_t p2 = 0; p2 < lattice.size(); ++p2)
    {
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));

		  for (pos_t p3 = ((bra_neq_ket) ? 0 : std::min(p1, p2)); p3 < lattice.size(); ++p3)
      {
		    if(p1 == p2 && p1 == p3)
          continue;

        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
        std::vector<std::vector<pos_t> > num_labels;

        for (pos_t p4 = 0; p4 < lattice.size(); ++p4)
        {
          if(std::max(p1,p2)  > std::max(p3,p4)  && std::max(p3,p4) > p2)
            continue;
          if(std::max(p1,p2) == std::max(p3,p4) && p1 < p4)
            continue;
		      if(p3 > p1 && p2 > p4)
            continue;
		      if(p1 == p4 && p3 == p2)
            continue;

          std::vector<pos_t> positions = {p1, p3, p4, p2};
          std::vector<pos_t> positionsORD = {p1, p2, p3, p4};

          // Loop over operator terms that are measured synchronously and added together
          typename MPS<Matrix, SymmGroup>::scalar_type value = 0;
          bool measured = false;
          for (std::size_t synop = 0; synop < operator_terms.size(); ++synop) {
            tag_vec operators(4);
            operators[0] = operator_terms[synop].first[0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
            operators[1] = operator_terms[synop].first[1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];
            operators[2] = operator_terms[synop].first[2][lattice.get_prop<typename SymmGroup::subcharge>("type", p3)];
            operators[3] = operator_terms[synop].first[3][lattice.get_prop<typename SymmGroup::subcharge>("type", p4)];
            term_descriptor term = generate_mpo::arrange_operators(positions, operators, tag_handler_local);
            // check if term is allowed by symmetry
            //if(not measurements_details::checkpg<SymmGroup>()(term, tag_handler_local, lattice))
            //    continue;
            measured = true;
            MPO<Matrix, SymmGroup> mpo = generate_mpo::sign_and_fill(term, identities, fillings, tag_handler_local, lattice);
            //value += operator_terms[synop].second * expval(bra_mps, ket_mps, mpo);
            value += (this->cast_to_real) ? maquis::real(operator_terms[synop].second * expval(bra_mps, ket_mps_local, mpo)) 
                                          : operator_terms[synop].second * expval(bra_mps, ket_mps_local, mpo);
          }

          if(measured)
          {
              dct.push_back(value);
              num_labels.push_back(order_labels(lattice, positionsORD));
          }
        }

        std::vector<std::string> lbt = label_strings(num_labels);

        // save results and labels
        #ifdef MAQUIS_OPENMP
        #pragma omp critical
        #endif
        {
          this->vector_results.reserve(this->vector_results.size() + dct.size());
          std::copy(dct.rbegin(), dct.rend(), std::back_inserter(this->vector_results));
          this->labels.reserve(this->labels.size() + dct.size());
          std::copy(lbt.rbegin(), lbt.rend(), std::back_inserter(this->labels));
          this->labels_num.reserve(this->labels_num.size() + dct.size());
          std::copy(num_labels.rbegin(), num_labels.rend(), std::back_inserter(this->labels_num));
        }
      }
    }
  }
  
private:
  Lattice lattice;
  std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
  positions_type positions_first;
  tag_vec identities, fillings;
  std::vector<scaled_bond_term> operator_terms;
  std::string bra_ckp;
};

} // namespace measurements

#endif // MEASUREMENTS_TAGGED_NRANKRDM_REL_H