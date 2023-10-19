/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_TAGGED_NRANKRDM_SU2U1_H
#define MEASUREMENTS_TAGGED_NRANKRDM_SU2U1_H

#include <algorithm>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>
#include "dmrg/block_matrix/symmetry/nu1pg.h"
#include "dmrg/models/measurement.h"
#include "dmrg/utils/checks.h"
#include "dmrg/models/chem/su2u1/term_maker.h"
#include "dmrg/models/chem/transform_symmetry.hpp"
#include "measurements_details.h"

namespace measurements {

template <class Matrix, class SymmGroup>
class TaggedNRankRDM<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup> > 
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
  typedef TermMakerSU2<Matrix, SymmGroup> TM;

public:
  TaggedNRankRDM(std::string const& name_, const Lattice & lat, std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                 typename TM::OperatorCollection const & op_collection_, positions_type const& positions_ = positions_type(),
                 std::string const& ckp_ = std::string(""))
    : base(name_), lattice(lat), tag_handler(tag_handler_), op_collection(op_collection_), positions_first(positions_),
      identities(op_collection.ident.no_couple), fillings(op_collection.fill.no_couple), bra_ckp(ckp_)
  {
    pos_t extent = lattice.size();
    if (positions_first.size() == 0)
    {
        positions_first.resize(extent);
        std::iota(positions_first.begin(), positions_first.end(), 0);
    }
    //this->cast_to_real = is_hermitian_meas(ops[0]);
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

    maquis::cout << " measuring in su2 version of tagged_nrank " << std::endl;

    if (this->name() == "oneptdm" || this->name() == "transition_oneptdm")
        measure_correlation(bra_mps, ket_mps);
    else if (this->name() == "twoptdm" || this->name() == "transition_twoptdm")
        measure_2rdm(bra_mps, ket_mps);
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
    #pragma omp parallel for schedule(dynamic) firstprivate(bra_mps, ket_mps_local)
    #endif
    for (std::size_t i = 0; i < positions_first.size(); ++i) {
      pos_t p1 = positions_first[i];
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));
      std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
      std::vector<std::vector<pos_t> > num_labels;
      for (pos_t p2 = (bra_neq_ket ? 0 : p1); p2 < lattice.size(); ++p2)
      {
        std::vector<pos_t> positions = {p1, p2};
        std::vector<term_descriptor> terms;
        if (p1 != p2) {
            // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
            terms.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, op_collection.ident.no_couple, std::sqrt(2.), p1, p2, op_collection.create.couple_down, op_collection.create.fill_couple_up,
                                                  op_collection.destroy.couple_down, op_collection.destroy.fill_couple_up, lattice
            ));
        }
        else {
            term_descriptor term;
            term.coeff = 1.;
            term.push_back( std::make_pair(p1, op_collection.count.no_couple[lattice.get_prop<typename SymmGroup::subcharge>("type", p1)]) );
            terms.push_back(term);
        }

        // check if term is allowed by symmetry
        if(not measurements_details::checkpg<SymmGroup>()(terms[0], tag_handler_local, lattice))
               continue;

        generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpo_m(lattice, op_collection.ident.no_couple, op_collection.ident_full.no_couple,
                                                              op_collection.fill.no_couple, tag_handler_local, terms);
        MPO<Matrix, SymmGroup> mpo = mpo_m.create_mpo();
        typename MPS<Matrix, SymmGroup>::scalar_type value = expval(bra_mps, ket_mps_local, mpo);

        dct.push_back(value);

        // use lattice.get_prop to save the positions already in correct order
        // thus label_strings below does not need to reorder the indices anymore
        num_labels.push_back(order_labels(lattice, positions));
      }

      // no need to reorder the labels here anymore
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

    // get all the labels ahead of the measurement and initialise the result arrays with the correct size
    auto indices = measurements_details::iterate_nrdm<2>(lattice.size(), bra_neq_ket);
    this->labels_num.resize(indices.size());
    this->labels.resize(indices.size());
    this->vector_results.resize(indices.size());

    #ifdef MAQUIS_OPENMP
    #pragma omp parallel for schedule(dynamic) firstprivate(ket_mps_local, bra_mps)
    #endif
    for (int i = 0; i < indices.size(); i++)
    {
      auto&& positions = indices[i];
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));
      std::vector<term_descriptor> terms = SpinSumSU2<Matrix, SymmGroup>::V_term(1., positions[0], positions[1], positions[2], positions[3], op_collection, lattice);
      // save labels
      auto&& num_labels = order_labels(lattice, positions);
      std::string lbt = label_string(num_labels);
      this->labels[i] = lbt;
      this->labels_num[i] = num_labels;

      // check if term is allowed by symmetry
      if(not measurements_details::checkpg<SymmGroup>()(terms[0], tag_handler_local, lattice)) {
          this->vector_results[i] = 0.;
          continue;
      }

      generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpo_m(lattice, op_collection.ident.no_couple, op_collection.ident_full.no_couple,
                                                              op_collection.fill.no_couple, tag_handler_local, terms);
      MPO<Matrix, SymmGroup> mpo = mpo_m.create_mpo();
      typename MPS<Matrix, SymmGroup>::scalar_type value = expval(bra_mps, ket_mps_local, mpo);

      // save results
      this->vector_results[i] = value;
    }
  }

private:
  Lattice lattice;
  std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
  typename TM::OperatorCollection op_collection;
  positions_type positions_first;
  tag_vec identities, fillings;
  std::string bra_ckp;
};

} // namespace measurements

#endif // MEASUREMENTS_TAGGED_NRANKRDM_SU2U1_H
