/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_TAGGED_NRANKRDM_TWOU1_H
#define MEASUREMENTS_TAGGED_NRANKRDM_TWOU1_H

namespace measurements {

template <class Matrix, class SymmGroup, class = void>
class TaggedNRankRDM : public measurement<Matrix, SymmGroup> {
  // Types declaration
  typedef measurement<Matrix, SymmGroup> base;
  typedef typename Model<Matrix, SymmGroup>::term_descriptor term_descriptor;
  typedef Lattice::pos_t pos_t;
  typedef std::vector<pos_t> positions_type;
  typedef typename base::op_t op_t;
  typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
  typedef typename base::value_type value_type;
  typedef std::vector<tag_type> tag_vec;
  typedef std::vector<tag_vec> bond_term;
  typedef std::pair<std::vector<tag_vec>, value_type> scaled_bond_term;

public:
  /** @brief Class constructor */
  TaggedNRankRDM(std::string const& name_, const Lattice & lat, std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
                 tag_vec const & identities_, tag_vec const & fillings_, std::vector<scaled_bond_term> const& ops_,
                 bool half_only_, positions_type const& positions_ = positions_type(), std::string const& ckp_ = std::string(""))
    : base(name_), lattice(lat), tag_handler(tag_handler_), positions_first(positions_), identities(identities_),
      fillings(fillings_), operator_terms(ops_), half_only(half_only_), bra_ckp(ckp_)
  {
    pos_t extent = lattice.size();
    // the default setting is only required for "measure_correlation"
    if (positions_first.size() == 0 && operator_terms[0].first.size() == 2)
    {
        positions_first.resize(extent);
        std::iota(positions_first.begin(), positions_first.end(), 0);
    }
    //this->cast_to_real = is_hermitian_meas(ops[0]);
    this->cast_to_real = false;
  }

  /** @brief Calculates the RDM elements */
  void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
  {
    // Preliminary operations
    this->vector_results.clear();
    this->labels.clear();
    MPS<Matrix, SymmGroup> bra_mps;
    if (bra_ckp != "")
    {
      if(boost::filesystem::exists(bra_ckp))
      {
        // Do symmetry check on the bra checkpoint and eventually transform
        // check point group
        // or boost::is_same<HasPG<SymmGroup>, std::true_type>::value?
        if (maquis::checks::has_pg(bra_ckp) != symm_traits::HasPG<SymmGroup>::value)
            throw std::runtime_error("Bra checkpoint " + bra_ckp + "has the wrong point group symmetry.");
        // Check if the bra checkpoint has SU2 symmetry, if so, transform it
        std::regex su2_regex("^su2u1");
        std::string bra_sym = maquis::checks::detail::get_symmetry(bra_ckp);
        if (std::regex_search(bra_sym,su2_regex))
#if (defined(HAVE_SU2U1) || defined(HAVE_SU2U1PG))
        {
          typedef typename boost::mpl::if_<symm_traits::HasPG<SymmGroup>, SU2U1PG, SU2U1>::type SU2Symm;
          MPS<Matrix, SU2Symm> su2_mps;
          load(bra_ckp, su2_mps);
          int N = SU2Symm::particleNumber(su2_mps[su2_mps.size()-1].col_dim()[0].first);
          int TwoS = SU2Symm::spin(su2_mps[su2_mps.size()-1].col_dim()[0].first);
          int Nup = (N + TwoS) / 2;
          int Ndown = (N - TwoS) / 2;
          BaseParameters parms_tmp = chem::detail::set_2u1_parameters(su2_mps.size(), Nup, Ndown);
          //parms_tmp.set("MEASURE[ChemEntropy]", 1);
          Model<Matrix, SymmGroup> model_tmp(lattice, parms_tmp);
          bra_mps = transform_mps<Matrix, SU2Symm>()(su2_mps, Nup, Ndown);
        }
#else
          throw std::runtime_error("SU2U1 support has not been compiled, cannot transform the MPS from SU2U1 group.");
#endif
      else
        load(bra_ckp, bra_mps);
      }
      else
        throw std::runtime_error("The bra checkpoint file " + bra_ckp + " was not found\n");
    }
    maquis::cout << " measuring in 2u1 version of tagged_nrank " << std::endl;
    //
    if (operator_terms[0].first.size() == 2)
        measure_correlation(bra_mps, ket_mps);
    else if (operator_terms[0].first.size() == 4)
        measure_nrdm<2>(bra_mps, ket_mps);
    else if (operator_terms[0].first.size() == 6)
        measure_nrdm<3>(bra_mps, ket_mps);
    else if (operator_terms[0].first.size() == 8)
        measure_nrdm<4>(bra_mps, ket_mps);
    else
        throw std::runtime_error("correlation measurements at the moment supported with 2, 4, 6 and 8 operators, size is "
                                  + boost::lexical_cast<std::string>(operator_terms[0].first.size()));
  }

protected:

  /** @brief Cloning method */
  measurement<Matrix, SymmGroup>* do_clone() const { return new TaggedNRankRDM(*this); }

  /** @brief Off-diagonal measurement (i.e., bra != ket) */
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
      for (pos_t p2 = bra_neq_ket ? 0 : p1; p2 < lattice.size(); ++p2)
      {
        std::vector<pos_t> positions{p1, p2};
        typename MPS<Matrix, SymmGroup>::scalar_type value = 0.0;
        for (std::size_t synop = 0; synop < operator_terms.size(); ++synop)
        {
            tag_vec operators(2);
            operators[0] = operator_terms[synop].first[0][lattice.get_prop<typename SymmGroup::subcharge>("type", p1)];
            operators[1] = operator_terms[synop].first[1][lattice.get_prop<typename SymmGroup::subcharge>("type", p2)];
            // check if term is allowed by symmetry
            term_descriptor term = generate_mpo::arrange_operators(positions, operators, tag_handler_local);
            if(measurements_details::checkpg<SymmGroup>()(term, tag_handler_local, lattice))
            {
                MPO<Matrix, SymmGroup> mpo = generate_mpo::sign_and_fill(term, identities, fillings, tag_handler_local, lattice);
                value += operator_terms[synop].second * expval(bra_mps, ket_mps_local, mpo);
            }
        }
        dct.push_back(value);
        num_labels.push_back(order_labels(lattice, positions));
      }
      //
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
    
  // Generic function for measuring all RDMs
  template<int N>
  void measure_nrdm(const MPS<Matrix, SymmGroup> & dummy_bra_mps, const MPS<Matrix, SymmGroup> & ket_mps)
  {
      // Test if a separate bra state has been specified bool bra_neq_ket = (dummy_bra_mps.length() > 0);
      bool bra_neq_ket = (dummy_bra_mps.length() > 0);
      //MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
      MPS<Matrix, SymmGroup> bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;
      MPS<Matrix, SymmGroup> ket_mps_local = ket_mps;
      // Obtain the total number of RDM elements and the list of all indices (eventually for a given slice)
      auto indices = measurements_details::iterate_nrdm<N>(lattice.size(), bra_neq_ket, positions_first);
      maquis::cout << "Number of total " << N << "-RDM elements measured: " << indices.size() << std::endl;
      // Prepare result arrays
      resize_results(indices.size());
      // Loop over all indices
      #ifdef MAQUIS_OPENMP
      #pragma omp parallel for schedule(dynamic) firstprivate(bra_mps, ket_mps_local)
      #endif
      for (int i = 0; i < indices.size(); i++)
      {
          auto&& positions = indices[i];
          // Prepare labels
          auto&& num_labels = order_labels(lattice, positions);
          std::string lbt = label_string(num_labels);
          this->labels[i] = lbt;
          this->labels_num[i] = num_labels;
          // MPS<Matrix, SymmGroup> ket_mps_local = ket_mps; // enable if you get pairing issues
          // Make a local copy of tag_handler since it can be modified by the MPO creator
          std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_local(new TagHandler<Matrix, SymmGroup>(*tag_handler));
          // Setup MPO and calculate the expectation value for a given indices set
          this->vector_results[i] = nrdm_expval(N, bra_mps, ket_mps_local, positions, tag_handler_local);
      } // iterator loop
  }

private:
  // Class members
  Lattice lattice;
  std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
  positions_type positions_first;
  tag_vec identities, fillings;
  std::vector<scaled_bond_term> operator_terms;
  bool half_only;
  std::string bra_ckp;

  // Resize labels and results, used before the measurements
  inline void resize_results(int size)
  {
      this->labels_num.resize(size);
      this->labels.resize(size);
      this->vector_results.resize(size);
  }

  // Obtain an expectation value for <bra|op|ket> for given n-RDM order and positions
  inline value_type nrdm_expval(std::size_t n, const MPS<Matrix, SymmGroup> & bra_mps, const MPS<Matrix, SymmGroup> & ket_mps,
              const std::vector<int> & positions, const std::shared_ptr<TagHandler<Matrix, SymmGroup> > & tag_handler_local)
  {
      assert(operator_terms.size() > 0);
      auto opsize = operator_terms[0].first.size();
      assert(n*2 == opsize);
      value_type result = 0.;
      // spin combo loop
      for (std::size_t synop = 0; synop < operator_terms.size(); ++synop)
      {
          tag_vec operators(opsize);
          for (std::size_t op = 0; op < opsize; op++)
              operators[op] = operator_terms[synop].first[op][lattice.get_prop<typename SymmGroup::subcharge>("type", positions[op])];
          // check if term is allowed by symmetry
          term_descriptor term = generate_mpo::arrange_operators(positions, operators, tag_handler_local);
          if(!measurements_details::checkpg<SymmGroup>()(term, tag_handler_local, lattice))
              return 0.;
          MPO<Matrix, SymmGroup> mpo = generate_mpo::sign_and_fill(term, identities, fillings, tag_handler_local, lattice);
          result += operator_terms[synop].second * expval(bra_mps, ket_mps, mpo);
      }
      return result;
  }
};

} // namespace measurements 

#endif // MEASUREMENTS_TAGGED_NRANKRDM_TWOU1_H