/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2019- Institute for Theoretical Physics, ETH Zurich
 *               2019 by Anna Kelemen <akelemen@ethz.ch>
 *               2020- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef LINSOLVER_H
#define LINSOLVER_H

#include <complex>
// #include <Eigen/Core>
// #include <Eigen/Dense>
// #include <Eigen/IterativeLinearSolvers>
// #include <unsupported/Eigen/IterativeSolvers>
// #include <Eigen/Eigenvalues>
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
            std::shared_ptr<block_matrix<Matrix, SymmGroup>> precond)
    : sp_(sp), parms_(parms), rhsMPS_(rhsMPS), shift_(shift), precond_(precond) //, isFolded_(false)
  {
    //if (params["pI_folded"] == "yes")
    //    isFolded_ = true;
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
  std::pair<energy_type, MPSTensorType> res() {
    maquis::cout << " Starting the iterative solution to the linear system" << std::endl;
    std::cout << " Norm of the rhs MPSTensor = " << ietl::two_norm(rhsMPS_) << std::endl;
    auto initError = ietl::two_norm(applyOperator(currentSolution_)-rhsMPS_);
    std::cout << " Initial linear system error = " << initError << std::endl;
    for (int iCycle = 0; iCycle < numberOfMacroIterations_; iCycle++) {
      if (parms_["linsystem_solver"] == "GMRES")
        gmres();
      else if (parms_["linsystem_solver"] == "MINRES")
        minres();
      else
        throw std::runtime_error("[linsystem_solver] parameter not recognized");
    }
    auto finalError = ietl::two_norm(applyOperator(currentSolution_)-rhsMPS_);
    std::cout << " Final linear system error = " << finalError << std::endl;
    // == Finalization ==
    auto tmp2 = applyOperator(currentSolution_);
    // ietl::mult(sp, x, tmp2, 0, false);
    ietl::mult(*sp_, currentSolution_, tmp2);
    ScalarType tmp3 = ietl::dot(currentSolution_, tmp2) / ietl::dot(currentSolution_, currentSolution_);
    en = maquis::real(tmp3);
    std::pair<energy_type, MPSTensorType> r0 = std::make_pair(en, currentSolution_);
    return r0;
  };

  /** @brief Default class destructor */
  ~LinSolver() = default;

protected:

  /**
   * @brief Solve the local linear system with the GMRES algorithm.
   * The implementation is based on Saad's book on iterative methods.
   */
  void gmres() {
    // Printing
    printHeader();
    maquis::cout << " --------------------------------- " << std::endl;
    maquis::cout << " Iteration  | Rel. error estimate  " << std::endl;
    maquis::cout << " --------------------------------- " << std::endl;
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
      maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                   << residual[iter] << std::endl;
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
      vec_type result(iter, 0.);
      for (int i = 0; i < iter; i++) {
        ScalarType s = 0.;
        for (int j = i+1; j < iter; j++)
          s = s + R(i, j)*result[j];
        result[i] = (y[i] - s)/R(i, i);
      }
      for (int iFinal = 0; iFinal < iter; iFinal++)
        currentSolution_ += result[iFinal]*vecSpace[iFinal];
    }
    printEndl();
  }

  /** @brief Solve the linear system with MINRES (based on the PyKry python library). */
  void minres() {
    // Printing
    printHeader();
    maquis::cout << " --------------------------------- " << std::endl;
    maquis::cout << " Iteration  | Rel. error estimate  " << std::endl;
    maquis::cout << " --------------------------------- " << std::endl;
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
      maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                   << residual[iter] << std::endl;
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

  /** @brief Overloading of user-defined conjugate function */
  static inline double localConj(double a) { return a; }

  /** @brief Conjugate function for complex numbers */
  static inline std::complex<double> localConj(std::complex<double> a) { return std::conj(a); }

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
  static void printEndl() {
    maquis::cout << " --------------------------------- " << std::endl;
    maquis::cout << std::endl;
    maquis::cout << std::fixed;
  }

  /** @brief Prints the header of the table */
  static void printHeader() {
    maquis::cout << std::endl;
    maquis::cout << " +-----------------------------------------+ " << std::endl;
    maquis::cout << " | ITERATIVE SOLUTION OF THE LINEAR SYSTEM | " << std::endl;
    maquis::cout << " +-----------------------------------------+ " << std::endl;
  }

  /** @brief Calculates the Givens rotation */
  static std::tuple<ScalarType, ScalarType, ScalarType> givens(ScalarType x0, ScalarType x1, ScalarType in0, ScalarType in1) {
    double cos;
    ScalarType sin;
    std::tie(cos, sin) = calculateRotation(x0, x1);
    auto r = cos*x0 + sin*x1;
    auto out0 = cos*in0 + sin*in1;
    auto out1 = -localConj(sin)*in0+cos*in1;
    return std::make_tuple(out0, out1, r);
  }


  /* Private members */
  std::shared_ptr<SiteProblem<Matrix, SymmGroup>> sp_;       // Pointer to the site problem representing the linear system.
  std::shared_ptr<block_matrix<Matrix, SymmGroup>> precond_; // Pointer to the preconditioner.
  const MPSTensorType& rhsMPS_;                              // Reference to the MPS representing the RHS of the local linear system.
  BaseParameters& parms_;                                    // Parameter container.
  ScalarType shift_;                                         // Shift to apply to the Hamiltonian.
  energy_type en;                                            // CHECK IF NEEDED
  MPSTensorType currentSolution_;                            // Stores the current approximation to the solution of the linear system.
  int numberOfMacroIterations_;                              // Number of restarts for the solution of the linear system.
  int krylovDim_;                                            // Maximum dimension of the Krylov space.
  RealType gmresTol_;                                        // Convergence threshold for the iterative solution to the linear system.
  RealType rhsNorm_;                                         // Norm of the rhs term.
  static constexpr double zeroThresh_ = 1.0E-16;             // Numerical zero
  // bool isFolded_;
};

#endif
