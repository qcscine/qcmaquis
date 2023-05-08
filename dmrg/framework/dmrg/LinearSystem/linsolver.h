/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef LINSOLVER_H
#define LINSOLVER_H

#include <complex>
#include <tuple>
// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/IterativeSolvers>
// #include <Eigen/Eigenvalues>
#include <boost/numeric/bindings/lapack.hpp>
#include "linsolver_helper.h"
#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/siteproblem.h"

template<class Matrix, class SymmGroup>
class LinSolver
{
  using energy_type = typename MPSTensor<Matrix, SymmGroup>::magnitude_type;
  using ScalarType = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using RealType = typename MPSTensor<Matrix, SymmGroup>::real_type;
  // using mat_type = typename Eigen::Matrix<ScalarType, Eigen::Dynamic, Eigen::Dynamic>;
  // using vec_type = Eigen::Matrix<ScalarType, Eigen::Dynamic, 1>;
  using mat_type = Matrix;
  using vec_type = alps::numeric::vector<ScalarType>;
  using MPSTensorType = MPSTensor<Matrix, SymmGroup>;
  using basis_type = std::vector<MPSTensorType>;
  using GivensType = typename LinSolverHelper::Givens<ScalarType>;
  using BlockMatrixType = block_matrix<Matrix, SymmGroup>;

public:
  /**
   * @brief Constructor
   * @param sp Pointer to the SiteProblem representing the lhs of the linear system.
   * @param initialMPS Initial MPS for the site on which the sweep is centered.
   * @param rhsMPS MPS associated with the rhs term.
   * @param shift shift for the solution of the linear system.
   * @param parms parameter container.
   * @param precond pointer to the preconditioner. If set, applies the preconditioning.
   */
  LinSolver(std::shared_ptr<SiteProblem<Matrix, SymmGroup>> sp, const MPSTensorType& initialMPS,
            const MPSTensorType& rhsMPS, ScalarType shift, BaseParameters & parms,
            std::shared_ptr<block_matrix<Matrix, SymmGroup>> precond, bool verbose)
    : sp_(sp), parms_(parms), rhsMPS_(rhsMPS), shift_(shift), precond_(precond), verbose_(verbose)
    //, isFolded_(false)
  {
    if (parms_["linsystem_init"] == "zero")
      currentSolution_ = 0.*initialMPS;
    else
      currentSolution_ = initialMPS;
    //rhs.conjugate_inplace();
    // Parameters that are specific of the solution of the linear system.
    numberOfMacroIterations_ = parms_["linsystem_max_it"].as<int>();
    gmresTol_ = parms_["linsystem_tol"];
    krylovDim_ = parms_["linsystem_krylov_dim"];
    rhsNorm_ = ietl::two_norm(rhsMPS_);
  }

  /** @brief Solves the linear system */
  std::tuple<energy_type, energy_type, MPSTensorType> res() {
    int prec = maquis::cout.precision();
    if (verbose_) {
      maquis::cout.precision(15);
      maquis::cout << std::endl;
      auto tmp2 = applyOperator(currentSolution_);
      auto initError = ietl::two_norm(tmp2-rhsMPS_);
      maquis::cout << " Initial ||Ax - b|| norm =  " << initError << std::endl;
      maquis::cout << std::endl;
    }
    for (int iCycle = 0; iCycle < numberOfMacroIterations_; iCycle++) {
      if (parms_["linsystem_solver"] == "GMRES")
        gmres();
      else if (parms_["linsystem_solver"] == "MINRES")
        minres();
      else
        throw std::runtime_error("[linsystem_solver] parameter not recognized");
    }
    auto tmp3 = applyOperator(currentSolution_);
    auto finalError = ietl::two_norm(tmp3-rhsMPS_);
    if (verbose_) {
      maquis::cout << std::endl;
      maquis::cout << " Final ||Ax - b|| norm =  " << finalError << std::endl;
    }
    // == Finalization ==
    // ietl::mult(sp, x, tmp2, 0, false);
    ietl::mult(*sp_, currentSolution_, tmp3);
    auto en = maquis::real(ietl::dot(currentSolution_, tmp3) / ietl::dot(currentSolution_, currentSolution_));
    if (verbose_) {
      maquis::cout << " Final energy = " << en << std::endl;
      maquis::cout << std::endl;
    }
    maquis::cout.precision(prec);
    return std::make_tuple(en, finalError, currentSolution_);
  };

  /** @brief Default class destructor */
  ~LinSolver() = default;

protected:

  /**
   * @brief Solve the local linear system with the GMRES algorithm.
   * The implementation is based on Saad's book on iterative methods.
   */
  void gmres() {
    if (verbose_) {
      maquis::cout << " ------------------------------------- " << std::endl;
      maquis::cout << " Iteration  | Rel. error estimate      " << std::endl;
      maquis::cout << " ------------------------------------- " << std::endl;
    }
    // Sets up the initial value of all parameters.
    int iter = 0;
    bool exit = false;
    std::vector<MPSTensorType> vecSpace, precondVecSpace;
    std::vector<double> residual;
    residual.reserve(krylovDim_);
    vecSpace.reserve(krylovDim_);
    if (precond_)
      precondVecSpace.reserve(krylovDim_);
    // vec_type y = vec_type::Zero(krylov_dim+1);
    vec_type y = vec_type(krylovDim_+1, 0.);
    std::vector<GivensType> givensRotations;
    givensRotations.reserve(krylovDim_);
    // mat_type H = mat_type::Zero(krylov_dim+1, krylov_dim);
    // mat_type R = mat_type::Zero(krylov_dim+1, krylov_dim);
    mat_type H(krylovDim_+1, krylovDim_, 0.);
    mat_type R(krylovDim_+1, krylovDim_, 0.);
    MPSTensorType initialError = rhsMPS_ - applyOperator(currentSolution_), preconditionedError;
    double initialErrorNorm;
    if (precond_) {
      preconditionedError = initialError;
      precond(preconditionedError);
      initialErrorNorm = std::sqrt(maquis::real(ietl::dot(initialError, preconditionedError)));
    }
    else {
      initialErrorNorm = ietl::two_norm(initialError);
    }
    residual.push_back(initialErrorNorm);
    initialError /= initialErrorNorm;
    if (residual[0] > zeroThresh_) {
      if (precond_) {
        preconditionedError /= initialErrorNorm;
        vecSpace.push_back(preconditionedError);
        precondVecSpace.push_back(initialError);
      }
      else {
        vecSpace.push_back(initialError);
      }
    }
    y[0] = residual[0];
    // == MAIN LOOP ==
    while (residual[iter] > gmresTol_ && iter < krylovDim_-1 && !exit) {
      if (verbose_) {
        maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                     << residual[iter] << std::endl;
      }
      // Begin of the Arnoldi part
      auto Av = applyOperator(vecSpace[iter]);
      if (iter > 0) {
          H(iter-1, iter) = H(iter, iter-1);
          if (precond_)
              Av -= H(iter, iter-1)*precondVecSpace[iter-1];
          else
              Av -= H(iter, iter-1)*vecSpace[iter-1];
      }
      auto alpha = ietl::dot(vecSpace[iter], Av);
      if (precond_)
          Av -= alpha*precondVecSpace[iter];
      else
          Av -= alpha*vecSpace[iter];
      H(iter, iter) += alpha;
      MPSTensorType pAv;
      if (precond_) {
          pAv = Av;
          precond(pAv);
          H(iter+1, iter) = std::sqrt(ietl::dot(Av, pAv));
      }
      else {
          H(iter+1, iter) = ietl::two_norm(Av);
      }
      // Check if it's zero norm
      if (std::abs(H(iter+1, iter)) < zeroThresh_) {
          exit = true;
      }
      else {
          if (precond_) {
              vecSpace.push_back(pAv/H(iter+1, iter));
              precondVecSpace.push_back(Av/H(iter+1, iter));
          }
          else {
              vecSpace.push_back(Av/H(iter+1, iter));
          }
      }
      // Solution of the linear system
      for (int iRow = 0; iRow < iter+2; iRow++)
          R(iRow, iter) = H(iRow, iter);
      for (int iPair = 0; iPair < iter; iPair++)
          std::tie(R(iPair, iter), R(iPair+1, iter)) = givensRotations[iPair].apply(R(iPair, iter), R(iPair+1, iter));
      givensRotations.emplace_back(GivensType(R(iter, iter), R(iter+1, iter)));
      std::tie(R(iter, iter), R(iter+1, iter)) = givensRotations[iter].apply(R(iter, iter), R(iter+1, iter));
      std::tie(y[iter], y[iter+1]) = givensRotations[iter].apply(y[iter], y[iter+1]);
      residual.push_back(std::abs(y[iter+1])/rhsNorm_);
      iter += 1;
    }
    // Final back-substitution
    if (iter != 0) {
      //vec_type result = R.block(0, 0, iter, iter).colPivHouseholderQr().solve(y.head(iter));
      //vec_type result = vec_type::Zero(iter);
      mat_type smallerMatrix(iter, iter, 0.);
      vec_type smallerVector(iter);
      for (int i = 0; i < iter; i++) {
        for (int j = 0; j < iter; j++)
          smallerMatrix(i, j) = R(i, j);
        smallerVector[i] = y[i];
      }
      /*
      vec_type result(iter, 0.);
      for (int i = 0; i < iter; i++) {
        ScalarType s = 0.;
        for (int j = i+1; j < iter; j++)
          s = s + R(i, j)*result[j];
        result[i] = (y[i] - s)/R(i, i);
      }
      */
      auto info = boost::numeric::bindings::lapack::gels(smallerMatrix, smallerVector);
      if (info != 0)
        throw std::runtime_error("Error in the solution of the linear systen");
      for (int iFinal = 0; iFinal < iter; iFinal++)
        currentSolution_ += smallerVector[iFinal]*vecSpace[iFinal];
        //currentSolution_ += result[iFinal]*vecSpace[iFinal];
      // DEBUG
      // vec_type diff(iter, 0.);
      // for (int iRow = 0; iRow < iter; iRow++) {
      //   for (int iCol = 0; iCol < iter; iCol++)
      //     diff(iRow) += R(iRow, iCol)*smallerVector[iCol];
      //   std::cout << diff(iRow) - y[iRow] << std::endl;
      // }
    }
    printEndl();
  }

  /** @brief Solve the linear system with MINRES (based on the PyKry python library). */
  void minres() {
    // Printing
    if (verbose_) {
      maquis::cout << " ------------------------------------- " << std::endl;
      maquis::cout << " Iteration  | Rel. error estimate      " << std::endl;
      maquis::cout << " ------------------------------------- " << std::endl;
    }
    // Sets up the initial value of all parameters.
    int iter = 0;
    bool exit = false;
    std::vector<MPSTensorType> vecSpace, precondVecSpace;
    std::vector<double> residual;
    residual.reserve(krylovDim_);
    vecSpace.reserve(krylovDim_);
    if (precond_)
      precondVecSpace.reserve(krylovDim_);
    // mat_type H = mat_type::Zero(krylov_dim+1, krylov_dim);
    mat_type H(krylovDim_+1, krylovDim_, 0.);
    auto v = applyOperator(rhsMPS_);
    residual.push_back(ietl::two_norm(rhsMPS_ - v));
    ScalarType vnorm = ietl::dot(rhsMPS_, v);
    if (std::norm(vnorm) > zeroThresh_) {
      precondVecSpace.push_back(rhsMPS_/vnorm);
      vecSpace.push_back(v/vnorm);
    }
    auto W = std::vector<MPSTensorType>(2);
    std::vector<ScalarType> y(2);
    y[0] = residual[0];
    auto yk = 0.*rhsMPS_;
    GivensType G1, G2;
    MPSTensorType MAv;
    // == MAIN LOOP ==
    while (residual[iter] > gmresTol_ && iter < krylovDim_-1 && !exit) {
      if (verbose_) {
        maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                     << residual[iter] << std::endl;
      }
      auto Av = applyOperator(vecSpace[iter]);
      if (iter > 0) {
        H(iter-1, iter) = H(iter, iter-1);
        if (precond_)
          Av -= H(iter, iter-1)*precondVecSpace[iter-1];
        else
          Av -= H(iter, iter-1)*vecSpace[iter-1];
      }
      auto alpha = ietl::dot(vecSpace[iter], Av);
      if (precond_)
        Av -= alpha*precondVecSpace[iter];
      else
        Av -= alpha*vecSpace[iter];
      H(iter, iter) = alpha;
      if (precond_) {
          MAv = Av;
          precond(MAv);
          H(iter+1, iter) = ietl::dot(Av, MAv);
      } else {
          H(iter+1, iter) = ietl::dot(Av, Av);
      }
      // Check if it's zero norm. If this is the case, exit.
      if (std::abs(H(iter+1, iter)) < zeroThresh_) {
        exit = true;
      }
      else {
        if (precond_) {
          precondVecSpace.push_back(Av/H(iter+1, iter));
          vecSpace.push_back(MAv/H(iter+1, iter));
        }
        else {
          vecSpace.push_back(Av/H(iter+1, iter));
        }
      }
      iter += 1;
      // vec_type R = vec_type::Zero(4);
      vec_type R(4, 0.);
      R(1) = maquis::real(H(iter-1, iter));
      if (G1.isActivated())
        std::tie(R(0), R(1)) = G1.apply(R(0), R(1));
      R(2) = maquis::real(H(iter, iter));
      R(3) = maquis::real(H(iter+1, iter));
      if (G2.isActivated())
        std::tie(R(1), R(2)) = G1.apply(R(1), R(2));
      G1 = G2;
      G2 = GivensType(R(2), R(3));
      R(2) = G2.getR();
      R(3) = 0.0;
      std::tie(y[0], y[1]) = G2.apply(y[0], y[1]);
      auto z = vecSpace[iter]/R(2);
      if (iter > 2)
        z -= R(0)*W[0]/R(2);
      if (iter > 1)
        z -= R(1)*W[1]/R(2);
      W[0] = W[1];
      W[1] = z;
      yk = yk + y[0]*z;
      y[0] = y[1];
      y[1] = 0.;
      residual.push_back(std::abs(y[0])/rhsNorm_);
    }
    currentSolution_ += yk;
    printEndl();
  }

public:
  /**
   * @brief Apply the shifted operator onto the MPS.
   * @param inputVec lhs vector
   * @return MPSTensorType rhs vector
   */
  MPSTensorType applyOperator(const MPSTensorType& inputVec) const {
    MPSTensorType ret; //, retSquared;
    ietl::mult(*sp_, inputVec, ret);
    //if (isFolded_) {
    //    ietl::mult(sp, inputVec, retSquared, 0, true);
    //    if (params["lin_alg"] == "feast")
    //        ret = retSquared - 2.*maquis::real(Z)*ret + inputVec*std::norm(Z);
    //    else
    //        ret = retSquared - 2.*sigma*ret + inputVec*sigma*sigma;
    //}
    //else {
    ret = ret - shift_*inputVec;
    return ret;
  }

private:

  /** @brief Preconditioner */
  void precond(MPSTensorType& inputVec) const {
    auto& data = inputVec.data();
    ScalarType denom;
    assert(shape_equal(data, *precond_));
    for (int b = 0; b < data.n_blocks(); b++) {
      for (size_t i = 0; i < num_rows(data[b]); ++i) {
        for (size_t j = 0; j < num_cols(data[b]); ++j) {
          denom = (precond_->operator[](b)(i, j) - shift_);
          if (std::fabs(denom) > 1.0E-10)
            data[b](i, j) /= std::fabs(denom);
        }
      }
    }
  }

  /** @brief Just prints a line for the table of the results */
  void printEndl() {
    if (verbose_) {
      maquis::cout << " ------------------------------------- " << std::endl;
      maquis::cout << std::endl;
      maquis::cout << std::fixed;
    }
  }

  /* Private members */
  std::shared_ptr<SiteProblem<Matrix, SymmGroup>> sp_;       // Pointer to the site problem representing the linear system.
  std::shared_ptr<block_matrix<Matrix, SymmGroup>> precond_; // Pointer to the preconditioner.
  const MPSTensorType& rhsMPS_;                              // Reference to the MPS representing the RHS of the local linear system.
  BaseParameters& parms_;                                    // Parameter container.
  ScalarType shift_;                                         // Shift to apply to the Hamiltonian.
  MPSTensorType currentSolution_;                            // Stores the current approximation to the solution of the linear system.
  int numberOfMacroIterations_;                              // Number of restarts for the solution of the linear system.
  int krylovDim_;                                            // Maximum dimension of the Krylov space.
  RealType gmresTol_;                                        // Convergence threshold for the iterative solution to the linear system.
  RealType rhsNorm_;                                         // Norm of the rhs term.
  static constexpr double zeroThresh_ = 1.0E-16;             // Numerical zero
  bool verbose_;                                             // Verbosity flag
};

#endif
