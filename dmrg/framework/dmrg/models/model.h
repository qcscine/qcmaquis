/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODELS_MODELS_H
#define MAQUIS_DMRG_MODELS_MODELS_H

#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/measurement.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/OperatorHandlers/TagHandler.h"
#include "dmrg/mp_tensors/mps_initializers.h"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/site_operator.h"

/// forward declaration
template <class Matrix, class SymmGroup>
class Measurements;

/// base type for all models
template <class Matrix, class SymmGroup>
class model_impl {
 public:
  using initializer_ptr = std::shared_ptr<mps_initializer<Matrix, SymmGroup>>;

  using table_type = TagHandler<Matrix, SymmGroup>;
  using table_ptr = std::shared_ptr<table_type>;
  using tag_type = typename table_type::tag_type;

  using term_descriptor = ::term_descriptor<typename Matrix::value_type>;
  using terms_type = typename std::vector<term_descriptor>;
  using op_t = typename operator_selector<Matrix, SymmGroup>::type;
  using measurements_type = boost::ptr_vector<measurement<Matrix, SymmGroup>>;
  using meas_with_results_type = std::pair<
      std::vector<std::vector<int>>, std::vector<typename Matrix::value_type>>;
  using results_map_type = std::map<std::string, meas_with_results_type>;

  using size_t = std::size_t;

  virtual ~model_impl() = default;

  virtual void update(BaseParameters const& p) = 0;

  virtual Index<SymmGroup> const& phys_dim(size_t type) const = 0;
  virtual op_t const& identity_matrix(size_t type) const {
    return operators_table()->get_op(identity_matrix_tag(type));
  }
  virtual tag_type identity_matrix_tag(size_t type) const = 0;
  virtual op_t const& filling_matrix(size_t type) const {
    return operators_table()->get_op(filling_matrix_tag(type));
  }
  virtual tag_type filling_matrix_tag(size_t type) const = 0;

  virtual typename SymmGroup::charge total_quantum_numbers(BaseParameters& parms
  ) const = 0;

  virtual terms_type const& hamiltonian_terms() const { return terms_; }
  virtual measurements_type measurements() const = 0;

  virtual op_t const& get_operator(std::string const& name, size_t type) const {
    return operators_table()->get_op(get_operator_tag(name, type));
  }
  virtual tag_type get_operator_tag(std::string const& name, size_t type)
      const = 0;

  virtual table_ptr operators_table() const = 0;

  virtual initializer_ptr initializer(Lattice const& lat, BaseParameters& parms)
      const;

  // optionally delay the assemly of the operator terms until the MPO is
  // actually created
  virtual void create_terms(){};

 protected:
  terms_type terms_;
};

/// model factory
template <class Matrix, class SymmGroup>
std::shared_ptr<model_impl<Matrix, SymmGroup>> model_factory(
    Lattice const& lattice, BaseParameters& parms
);

/// pimpl for Model
template <class Matrix, class SymmGroup>
class Model {
  using impl_type = model_impl<Matrix, SymmGroup>;
  using impl_ptr = std::shared_ptr<impl_type>;

 public:
  using initializer_ptr = typename impl_type::initializer_ptr;

  using table_type = typename impl_type::table_type;
  using table_ptr = typename impl_type::table_ptr;
  using tag_type = typename impl_type::tag_type;

  using term_descriptor = typename impl_type::term_descriptor;
  using terms_type = typename impl_type::terms_type;
  using op_t = typename impl_type::op_t;
  using measurements_type = typename impl_type::measurements_type;
  using meas_with_results_type = typename impl_type::meas_with_results_type;
  using results_map_type = typename impl_type::results_map_type;

  using size_t = typename impl_type::size_t;

  Model() = default;

  Model(Lattice const& lattice, BaseParameters& parms)
      : impl_(model_factory<Matrix, SymmGroup>(lattice, parms)) {}

  Model(impl_ptr impl) : impl_(impl) {}

  void update(BaseParameters const& p) { return impl_->update(p); }

  Index<SymmGroup> const& phys_dim(size_t type = 0) const {
    return impl_->phys_dim(type);
  }
  op_t const& identity_matrix(size_t type = 0) const {
    return impl_->identity_matrix(type);
  }
  tag_type identity_matrix_tag(size_t type = 0) const {
    return impl_->identity_matrix_tag(type);
  }
  op_t const& filling_matrix(size_t type = 0) const {
    return impl_->filling_matrix(type);
  }
  tag_type filling_matrix_tag(size_t type = 0) const {
    return impl_->filling_matrix_tag(type);
  }

  typename SymmGroup::charge total_quantum_numbers(BaseParameters& parms
  ) const {
    return impl_->total_quantum_numbers(parms);
  }

  terms_type const& hamiltonian_terms() const {
    return impl_->hamiltonian_terms();
  }
  measurements_type measurements() const { return impl_->measurements(); }

  op_t const& get_operator(std::string const& name, size_t type = 0) const {
    return impl_->get_operator(name, type);
  }
  tag_type get_operator_tag(std::string const& name, size_t type = 0) const {
    return impl_->get_operator_tag(name, type);
  }

  table_ptr operators_table() const { return impl_->operators_table(); }

  initializer_ptr initializer(Lattice const& lat, BaseParameters& parms) const {
    return impl_->initializer(lat, parms);
  }

  void create_terms() { impl_->create_terms(); }

 private:
  impl_ptr impl_;
};

#endif
