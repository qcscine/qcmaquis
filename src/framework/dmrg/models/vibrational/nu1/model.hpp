/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef MODELS_VIBRATIONAL_NU1_H
#define MODELS_VIBRATIONAL_NU1_H

#ifdef DMRG_VIBRATIONAL

#include <sstream>
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models//model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalIntegralParser.hpp"
#include "dmrg/models/model_helper.hpp"
#include "VibrationalRDM.h"

/**
 * @brief Class implementing the n-mode vibrational Hamiltonian
 *
 * This Hamiltonian relies on the SQ formalism introduced by Ove Christiansen in
 * JCP, 120, 2140 (2004) and, unlike the more common canonical quantization,
 * introduces a pair of SQ operators per modals *and* per mode. In this way,
 * it is possible to encode entirely generic potential operators that are
 * expressed in terms of entirely generic modal bases.
 */

template <class Matrix, int N>
class NMode : public model_impl<Matrix, NU1_template<N>> {
  // Types definition
  using NU1 = NU1_template<N>;
  using base = model_impl<Matrix, NU1>;
  using table_type = typename base::table_type;
  using table_ptr = typename base::table_ptr;
  using tag_type = typename base::tag_type;
  using term_descriptor = typename base::term_descriptor;
  using terms_type = typename std::vector<term_descriptor>;
  using op_t = typename base::op_t;
  using operators_type = typename std::vector<tag_type>;
  using measurements_type = typename base::measurements_type;
  using pos_t = typename Lattice::pos_t;
  using positions_type = typename std::vector<pos_t>;
  using value_type = typename Matrix::value_type;
  using charge_type = typename NU1::charge;

 public:
  /**
   * @brief Class constructor
   * @param lattice object representing the DMRG lattice
   * @param parameters container with the DMRG parameters
   * @param verbose if true, prints information regarding the Hamiltonian terms
   */
  NMode(const Lattice& lattice_, BaseParameters& parameters_, bool verbose)
      : lattice(lattice_),
        parameters(parameters_),
        tag_handler(new table_type()),
        phys_indexes(0) {
    // Loads in the relevant parameters
    num_modes = parameters["nmode_num_modes"];
    maxCouplingDegree = parameters["nmode_max_coupling"];
    lattice_size = parameters["L"];
    // == BUILDS ALL THE RELEVANT CHARGES ==
    // The ones we are interested in are:
    // - the empty charge
    // - the charges in which a given mode is populated by 1 quantum
    charge_type empty_state(0);
    std::vector<charge_type> excited_states(num_modes, 0);
    for (std::size_t idx = 0; idx < num_modes; idx++)
      excited_states[idx][idx] = 1;
    // Loads the physical indexes
    phys_indexes.reserve(num_modes);
    for (std::size_t idx = 0; idx < num_modes; idx++) {
      Index<NU1> phys;
      phys.insert(std::make_pair(empty_state, 1));
      phys.insert(std::make_pair(excited_states[idx], 1));
      phys_indexes.push_back(phys);
    }
    // Loads the vector with the site types
    siteTypes.reserve(lattice_size);
    for (int iSite = 0; iSite < lattice_size; iSite++)
      siteTypes.push_back(lattice.get_prop<int>("type", iSite));
    //
    // == DEFINITION OF THE ELEMENTARY OPERATORS ==
    // For each mode, we define the identity, the creation, the annihilation,
    // and the count operator. Note that we just need one operator per *mode*,
    // and not one per modal, since SQ operators acting on different modals will
    // anyway be mapped to different sites.
    std::vector<op_t> ident_op, create_op, destroy_op, count_op;
    ident.reserve(num_modes);
    create.reserve(num_modes);
    destroy.reserve(num_modes);
    count.reserve(num_modes);
    for (int idx = 0; idx < num_modes; idx++) {
      // Local operators
      op_t ident_op_loc, create_op_loc, destroy_op_loc, count_op_loc;
      ident_op_loc.insert_block(Matrix(1, 1, 1), empty_state, empty_state);
      ident_op_loc.insert_block(
          Matrix(1, 1, 1), excited_states[idx], excited_states[idx]
      );
      create_op_loc.insert_block(
          Matrix(1, 1, 1), empty_state, excited_states[idx]
      );
      destroy_op_loc.insert_block(
          Matrix(1, 1, 1), excited_states[idx], empty_state
      );
      count_op_loc.insert_block(
          Matrix(1, 1, 1), excited_states[idx], excited_states[idx]
      );
      // Updates the vectors
      ident_op.push_back(ident_op_loc);
      create_op.push_back(create_op_loc);
      count_op.push_back(count_op_loc);
      destroy_op.push_back(destroy_op_loc);
    }
    // Creates the final tags and update the table
    ident = modelHelper<Matrix, NU1>::register_all_types(
        ident_op, tag_detail::bosonic, tag_handler
    );
    create = modelHelper<Matrix, NU1>::register_all_types(
        create_op, tag_detail::bosonic, tag_handler
    );
    destroy = modelHelper<Matrix, NU1>::register_all_types(
        destroy_op, tag_detail::bosonic, tag_handler
    );
    count = modelHelper<Matrix, NU1>::register_all_types(
        count_op, tag_detail::bosonic, tag_handler
    );
    // Registers the hermitian pairs
    modelHelper<Matrix, NU1>::registerHermitianConjugates(
        create, destroy, tag_handler
    );
  }

  /** @brief Update the model with the new parameters */
  void update(BaseParameters const& p) {
    // TODO: update this->terms_ with the new parameters
    throw std::runtime_error("update() not yet implemented for this model.");
  }

  void create_terms() override {
    std::cout << "Parsing integral file" << std::endl;
    auto Hamiltonian_term =
        Vibrational::detail::NModeIntegralParser<double>(parameters, lattice);
    int hamiltonianSize = Hamiltonian_term.first.size();
    std::cout << "Processing Second-Quantization Hamiltonian" << std::endl;
    for (int iTerm = 0; iTerm < hamiltonianSize; iTerm++) {
      positions_type positions;
      operators_type operators;
      convertLineToOperators(
          Hamiltonian_term.first[iTerm], positions, operators
      );
      if (positions.size() / 2 <= maxCouplingDegree) {
        auto matrixElement =
            static_cast<value_type>(Hamiltonian_term.second[iTerm]);
        modelHelper<Matrix, NU1>::add_term(
            positions, operators, matrixElement, tag_handler, this->terms_
        );
      }
    }
    std::cout << "Second-Quantization Hamiltonian processed" << std::endl;
  }

  /** @brief Getter for the physical dimension of a given type */
  Index<NU1> const& phys_dim(size_t type) const { return phys_indexes[type]; }

  /** @brief Getter for the identity operator */
  tag_type identity_matrix_tag(size_t type) const { return ident[type]; }

  /** @brief Getter for the filling operator */
  tag_type filling_matrix_tag(size_t type) const {
    return identity_matrix_tag(type);
  }

  /**
   * @brief Gets the quantum number associated with the wfn
   *
   * The quantum number is generated by assigning 1 quantum to each mode.
   */
  typename NU1::charge total_quantum_numbers(BaseParameters& parms) const {
    return typename NU1::charge(std::vector<int>(num_modes, 1));
  }

  /**
   * @brief Gets the operator associated with a given string
   * @param name string describing the operator
   * @param type site type for which the operator is returned
   * @return tag_type tag associated with the requested operator
   */
  tag_type get_operator_tag(const std::string& name, size_t type) const {
    if (name == "n")
      return count[type];
    else if (name == "bdag")
      return create[type];
    else if (name == "b")
      return destroy[type];
    else if (name == "id")
      return ident[type];
    else if (name == "fill")
      return ident[type];
    else
      throw std::runtime_error("Operator not valid for this model.");
    return 0;
  }

  /** @brief Getter for the tag_handler */
  table_ptr operators_table() const { return tag_handler; }

  /**
   * @brief Measurement associated with the n-mode Hamiltonian class
   * Note that we currently support one- and two-modal RDM.
   */
  measurements_type measurements() const {
    measurements_type meas;
    // == One-modal RDM ==
    std::string name;
    std::vector<std::vector<operators_type>> ops;
    std::vector<float_t> coeffs;
    // One-mode RDM
    if (this->parameters.is_set("MEASURE[One Mode RDM]")) {
      name = "onemodeRDM";
      ops = {{create, destroy}};
      coeffs = {1.0};
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::AddAll
      ));
    }
    if (this->parameters.is_set("MEASURE[One Modal RDM]")) {
      // Calculate the one modal reduced density matrix for every mode
      // 1) operator 1 - a+a
      name = "onemodalRDM_00";
      ops = {{count}, {ident}};
      coeffs = {-1.0, 1.0};
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 2) operator a+a
      name = "onemodalRDM_11";
      ops = {{count}};
      coeffs = {1.0};
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
    }
    // == Two-modal RDM ==
    if (this->parameters.is_set("MEASURE[Two Modal RDM]")) {
      // 0-0 element of the two-particle RDM
      ops = {{ident, ident}, {count, ident}, {ident, count}, {count, count}};
      coeffs = {1.0, -1.0, -1.0, 1.0};
      name = "twomodeRDM_00";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 1-1 element of the two-particle RDM
      ops = {{count, ident}, {count, count}};
      coeffs = {1.0, -1.0};
      name = "twomodeRDM_11";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 1-2 element of the two-particle RDM
      ops = {{create, destroy}};
      coeffs = {1.0};
      name = "twomodeRDM_12";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 2-1 element of the two-particle RDM
      ops = {{destroy, create}};
      coeffs = {1.0};
      name = "twomodeRDM_21";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 2-2 element of the two-particle RDM
      ops = {{ident, count}, {count, count}};
      coeffs = {1.0, -1.0};
      name = "twomodeRDM_22";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
      // 3-3 element of the two-particle RDM
      ops = {{count, count}};
      coeffs = {1.0};
      name = "twomodeRDM_33";
      meas.push_back(new measurements::VibrationalRDM<Matrix, N>(
          lattice, name, ops, coeffs, this->tag_handler, this->ident,
          measurements::VibrationalRDMModality::ExcludeSame
      ));
    }
    return meas;
  }

 private:
  /**
   * @brief Converter of an input line to operators
   *
   * Note that, so far, we assume the lattice to be sorted so that we first have
   * all the modals associated with the first mode, then all the modals
   * assicated with the second mode and so on. Generic sortings are NYI.
   *
   * @param ham_term term to be parsed (array of integer numbers)
   * @param pos (output) vector with the position where operator are acting
   * @param ops (output) vector with the operators
   */
  template <class IntegralContainer>
  void convertLineToOperators(
      const IntegralContainer& ham_term, positions_type& pos,
      operators_type& ops
  ) {
    assert(ham_term.size() % 2 == 0);
    int jCont = 0;
    ops.reserve(ham_term.size());
    pos.reserve(ham_term.size());
    do {
      // Retrieves matrix element
      auto index = lattice.get_prop<int>(
          "absolutePositionInLattice", ham_term[2 * jCont] - 1,
          ham_term[2 * jCont + 1]
      );
      assert(index < lattice_size);
      pos.push_back(index);
      if (jCont % 2 == 0)  // since operators always come pairwise per mode
        ops.push_back(create[siteTypes[index]]);
      else
        ops.push_back(destroy[siteTypes[index]]);
      jCont += 1;
    } while (2 * jCont < ham_term.size() && ham_term[2 * jCont] != -1);
  }

 private:
  const Lattice& lattice;
  int lattice_size, num_modes, maxCouplingDegree;
  BaseParameters& parameters;
  std::vector<Index<NU1>> phys_indexes;
  std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler;
  operators_type ident, create, destroy, count;
  std::vector<int> siteTypes;
};

#endif  // DMRG_VIBRATIONAL

#endif
