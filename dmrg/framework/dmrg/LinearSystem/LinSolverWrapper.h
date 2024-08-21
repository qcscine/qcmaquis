/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#ifndef LINSOLVER_WRAPPER
#define LINSOLVER_WRAPPER

#include <Eigen/Core>
#include <Eigen/Sparse>

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/siteproblem.h"

/** @brief Methods that are implemented in Eigen for solving a linear system */
enum class EigenSolverType { BiCGSTAB, GMRES };

/** @brief Different types of preconditioners */
enum class PreconditionerType {
  IdentityPreconditioner,
  DiagonalPreconditioner
};

template <
    class MatrixWrapper, EigenSolverType EigenSolver,
    PreconditionerType Preconditioner>
class EigenSolverTraitClass;

template <class MatrixWrapper>
class EigenSolverTraitClass<
    MatrixWrapper, EigenSolverType::BiCGSTAB,
    PreconditionerType::IdentityPreconditioner> {
 public:
  using SolverType =
      Eigen::BiCGSTAB<MatrixWrapper, Eigen::IdentityPreconditioner>;
};

template <class MatrixWrapper>
class EigenSolverTraitClass<
    MatrixWrapper, EigenSolverType::BiCGSTAB,
    PreconditionerType::DiagonalPreconditioner> {
  using ScalarType = typename MatrixWrapper::Scalar;

 public:
  using SolverType =
      Eigen::BiCGSTAB<MatrixWrapper, Eigen::DiagonalPreconditioner<ScalarType>>;
};

template <class MatrixWrapper>
class EigenSolverTraitClass<
    MatrixWrapper, EigenSolverType::GMRES,
    PreconditionerType::IdentityPreconditioner> {
 public:
  using SolverType = Eigen::GMRES<MatrixWrapper, Eigen::IdentityPreconditioner>;
};

template <class MatrixWrapper>
class EigenSolverTraitClass<
    MatrixWrapper, EigenSolverType::GMRES,
    PreconditionerType::DiagonalPreconditioner> {
  using ScalarType = typename MatrixWrapper::Scalar;

 public:
  using SolverType =
      Eigen::GMRES<MatrixWrapper, Eigen::DiagonalPreconditioner<ScalarType>>;
};

/** Forward declaration of the SRIterativeMatrixWrapper class */
template <class Matrix, class SymmGroup>
class LinearSolverWrapper;

namespace Eigen {
namespace internal {
// Inherit traits from the most similar type in Eigen
template <class Matrix, class SymmGroup>
struct traits<LinearSolverWrapper<Matrix, SymmGroup>>
    : public traits<Eigen::Matrix<
          typename Matrix::value_type, Eigen::Dynamic, Eigen::Dynamic>> {};
}  // namespace internal
}  // namespace Eigen

/**
 * @class SRIterativeMatrixWrapper @file SRIterativeMatrixWrapper.h
 * @brief Class representing a matrix. It is used in iterative
 * @tparam ScalarType
 */
template <class Matrix, class SymmGroup>
class LinearSolverWrapper
    : public Eigen::EigenBase<LinearSolverWrapper<Matrix, SymmGroup>> {
 public:
  using Scalar = typename Matrix::value_type;
  using RealScalar = double;
  using StorageIndex = int;
  using Base = Eigen::EigenBase<LinearSolverWrapper<Matrix, SymmGroup>>;
  using SiteProblemType = SiteProblem<Matrix, SymmGroup>;
  using VectorType = MPSTensor<Matrix, SymmGroup>;
  using BlockMatrixType = block_matrix<Matrix, SymmGroup>;

  enum {
    ColsAtCompileTime = Eigen::Dynamic,
    MaxColsAtCompileTime = Eigen::Dynamic,
    IsRowMajor = false
  };

  /**
   * @brief Constructor taking the wave function gradients and their expectation
   * value.
   * @param gradients A matrix containing the determinant's gradient as columns
   * and parameters for rows.
   * @param gradientsExpectation The expectation value of the gradients over
   * determinants.
   * @param weights The weights of the determinants in the wave function.
   * @param mcmcSize The normalization factor of the wave function.
   */
  LinearSolverWrapper(
      std::shared_ptr<SiteProblemType> siteProblem, const VectorType& mpsTensor,
      std::shared_ptr<BlockMatrixType> preconditioner, Scalar z
  )
      : matrixFree(siteProblem), mpsReference(mpsTensor), shift(z) {
    if (preconditioner) populatePreconditioner(preconditioner);
  };

  /** @brief Getter for the rows of the wrapped matrix. */
  Eigen::Index rows() const { return mpsReference.num_elements(); }

  /** @brief Getter for the columns of the wrapped matrix. */
  Eigen::Index outerSize() const { return mpsReference.num_elements(); }

  /** @brief Updates the preconditioner */
  void populatePreconditioner(std::shared_ptr<BlockMatrixType> preconditioner) {
    const auto& data = mpsReference.data();
    assert(shape_equal(data, *preconditioner));
    preconditionValues.reserve(mpsReference.num_elements());
    for (int b = 0; b < data.n_blocks(); b++) {
      for (size_t i = 0; i < num_rows(data[b]); ++i) {
        for (size_t j = 0; j < num_cols(data[b]); ++j) {
          preconditionValues.push_back(preconditioner->operator[](b)(i, j));
        }
      }
    }
  }

  /** @brief Getter for the preconditioner coefficients */
  auto getPrecondValue(int idx) const {
    assert(idx < preconditionValues.size());
    return this->isPrecond() ? preconditionValues[idx]
                             : static_cast<Scalar>(1.);
  }

  /** @brief Checks if the object has a preconditioner */
  bool isPrecond() const { return preconditionValues.size() != 0; }

  /** @brief Getter for the columns of the wrapped matrix. */
  Eigen::Index cols() const { return mpsReference.num_elements(); }

  /** @brief Gets a copy of the reference MPS */
  VectorType getCopyRefMPSTensor() const { return mpsReference; }

  /** @brief Gets a reference to the underlying SiteProblem object */
  auto getSiteProblem() const { return matrixFree; }

  /** @brief Gets the shift */
  auto getShift() const { return shift; }

  /**
   * @brief Iterator class, needed for the preconditioner.
   * This code snippet was taken from:
   * https://stackoverflow.com/questions/51846432/diagonalpreconditioner-wrapper-on-eigen-3-3-5-for-matrix-free-sparse-solver-cg
   */
  class InnerIterator {
   public:
    // Types definition
    using Index = int;
    using MatrixReplacement = LinearSolverWrapper<Matrix, SymmGroup>;

    InnerIterator(const MatrixReplacement& mat, Index row)
        : mat(mat), has_val(false), row(row), col(row) {}

    /** @brief True if it's sitting on the diagonal */
    operator bool() { return row == col; }

    /** @brief Increases the operator by one */
    InnerIterator& operator++() {
      col++;
      has_val = false;
      return *this;
    }

    /** @brief Gette*/
    Index index() { return col; }

    /** @brief Black magic */
    Scalar value() {
      // Cache value since this function is called twice
      if (!has_val) {
        stored_val = mat.getPrecondValue(row);  // your implementation here
        has_val = true;
      }
      return stored_val;
    }

   private:
    const MatrixReplacement& mat;
    bool has_val;
    Index row, col;
    Scalar stored_val;
  };

  /**
   * @brief Product operator to use in expression templates.
   * @tparam Rhs The type of the vector multiplying the matrix.
   * @param rhs Product of this matrix wrapper with a vector.
   * @return An expression template describing a product operation.
   */
  template <typename Rhs>
  Eigen::Product<
      LinearSolverWrapper<Matrix, SymmGroup>, Rhs, Eigen::AliasFreeProduct>
  operator*(const Eigen::EigenBase<Rhs>& rhs) const {
    return Eigen::Product<
        LinearSolverWrapper<Matrix, SymmGroup>, Rhs, Eigen::AliasFreeProduct>(
        *this, rhs.derived()
    );
  }

 protected:
  std::vector<Scalar> preconditionValues;
  std::shared_ptr<SiteProblemType> matrixFree;
  const VectorType& mpsReference;
  Scalar shift;
};

namespace Eigen {
namespace internal {
/**
 * @brief Struct implementing a matrix-free Matrix-Vector product.
 * Here is where the meat is.
 * The iterative formulas are implemented in the function scalAndAddTo().
 */
template <typename Rhs, class Matrix, class SymmGroup>
struct generic_product_impl<
    LinearSolverWrapper<Matrix, SymmGroup>, Rhs, DenseShape, DenseShape,
    GemvProduct>
    : generic_product_impl_base<
          LinearSolverWrapper<Matrix, SymmGroup>, Rhs,
          generic_product_impl<LinearSolverWrapper<Matrix, SymmGroup>, Rhs>> {
  using MPSType = MPSTensor<Matrix, SymmGroup>;
  using WrapperType = LinearSolverWrapper<Matrix, SymmGroup>;
  using Scalar = typename Product<WrapperType, Rhs>::Scalar;

  /**
   * @brief Implement the effect of the matrix-vector multiplication. This is
   * the formation of the simga vector.
   * @tparam Destination The type of the result, a vector.
   * @param dest The vector resulting from the matrix-vector multiplication.
   * @param lhs The matrix-free wrapper acting as the matrix in the product.
   * @param rhs The vector onto which the matrix acts.
   * @param alpha A scalar. It is not used in the actual implementation.
   */
  template <typename Destination>
  static void scaleAndAddTo(
      Destination& dest, const WrapperType& lhs, const Rhs& rhs,
      const Scalar& alpha
  ) {
    // This method should implement "dst += alpha * lhs * rhs" inplace,
    // however, for iterative solvers, alpha is always equal to 1, so let's not
    // bother about it.
    assert(alpha == static_cast<Scalar>(1) && "scaling is not implemented");
    EIGEN_ONLY_USED_FOR_DEBUG(alpha);
    assert(
        dest.size() == rhs.size() &&
        "Destination vector has not the same size as input one in iterative SR."
    );

    // Converts the input vector in the MPS format
    auto inputMPS = lhs.getCopyRefMPSTensor();
    inputMPS.fillWithEigenVector(rhs);
    MPSType outputMPS = lhs.getSiteProblem()->apply(inputMPS);
    outputMPS = outputMPS - lhs.getShift() * inputMPS;
    auto outputVector = outputMPS.getEigenRepresentation();
    for (int iParam = 0; iParam < dest.size(); iParam++)
      dest(iParam) += outputVector(iParam);
  }
};

}  // namespace internal
}  // namespace Eigen

#endif  // LINSOLVER_WRAPPER