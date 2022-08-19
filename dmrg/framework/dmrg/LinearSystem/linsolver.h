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
class LinSolver{
public:
    using energy_type = typename MPSTensor<Matrix, SymmGroup>::magnitude_type;
    using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
    using real_type = typename MPSTensor<Matrix, SymmGroup>::real_type;
    // using mat_type = typename Eigen::Matrix<scalar_type, Eigen::Dynamic, Eigen::Dynamic>;
    // using vec_type = Eigen::Matrix<scalar_type, Eigen::Dynamic, 1>;
    using mat_type = Matrix;
    using vec_type = alps::numeric::vector<scalar_type>;
    using mps_tensor = typename MPSTensor<Matrix, SymmGroup>::MPSTensor;
    using basis_type = std::vector<mps_tensor>;
    using GivensType = typename LinSolverHelper::Givens<scalar_type>;

    /** @brief Class constructor */
    LinSolver(SiteProblem<Matrix, SymmGroup> const & sp_, MPSTensor<Matrix, SymmGroup> const & mpst_,
              MPSTensor<Matrix, SymmGroup> const & rhs_, scalar_type & Z_, BaseParameters & params_,
              block_matrix<Matrix, SymmGroup> const & precond)
        : sp(sp_), params(params_), rhs(rhs_), Z(Z_), sigma(0.), 
          precond_(precond), doPrecond_(false), isFolded_(false)
    {
        if (params["linsystem_precond"] == "yes")
            doPrecond_ = true;
        //if (params["pI_folded"] == "yes")
        //    isFolded_ = true;
        if (params["linsystem_init"] == "zero")
            x = 0.*mpst_;
        else
            x = mpst_;
        if (params.is_set("shift_exc"))
            sigma = params["shift_exc"];
        //rhs.conjugate_inplace();
        std::cout << "Initial norm " << ietl::two_norm(x) << std::endl;
    }
    
    std::pair<energy_type, mps_tensor> res() {
        const mps_tensor mpst_rhs = rhs;
        int gmres_max_it = params["linsystem_max_it"];
        real_type gmres_tol = params["linsystem_tol"];
        int krylov_dim = params["linsystem_krylov_dim"];
        // == Iterative Solver ==
        std::cout << " Rhs norm = " << ietl::two_norm(mpst_rhs) << std::endl;
        auto initError = ietl::two_norm(applyOperator(x)-mpst_rhs);
        std::cout << " Initial error = " << initError << std::endl;
        for (int iCycle = 0; iCycle < gmres_max_it; iCycle++) {
            if (params["linsystem_solver"] == "GMRES")
                gmres(sp, mpst_rhs, krylov_dim, gmres_tol);
            else if (params["linsystem_solver"] == "MINRES")
                minres(sp, mpst_rhs, krylov_dim, gmres_tol);
            else
                throw std::runtime_error("[linsystem_solver] parameter not recognized"); 
        }
        auto finalError = ietl::two_norm(applyOperator(x)-mpst_rhs);
        std::cout << " Final error = " << finalError << std::endl;
        // == Finalization ==
        mps_tensor tmp2 = applyOperator(x);
        // ietl::mult(sp, x, tmp2, 0, false);
        ietl::mult(sp, x, tmp2);
        scalar_type tmp3 = ietl::dot(x, tmp2) / ietl::dot(x, x);
        en = maquis::real(tmp3);
        std::pair<energy_type, mps_tensor> r0 = std::make_pair(en, x);
        return r0;
    };

    ~LinSolver() = default;

protected:

    /**
     * @brief Solve the local linear system with the GMRES algorithm.
     *
     * The implementation is based on Saad's book on iterative methods.
     *
     * @param sp site problem (representing the Hamiltotnian.
     * @param rhs right-hand side of the linear system.
     * @param krylov_dim maximum size of the Krylov subspace.
     * @param tol convergence tolerance for the whole GMRES procedure
     */
    void gmres(const SiteProblem<Matrix, SymmGroup> & sp, const mps_tensor& rhs,
               int krylov_dim, real_type tol) {
        // Printing
        printHeader();
        maquis::cout << " --------------------------------- " << std::endl;
        maquis::cout << " Iteration  | Rel. error estimate  " << std::endl;
        maquis::cout << " --------------------------------- " << std::endl;
        // Sets up the initial value of all parameters.
        int iter = 0;
        bool exit = false;
        std::vector<mps_tensor> vecSpace, precondVecSpace;
        std::vector<double> residual;
        residual.reserve(krylov_dim);
        vecSpace.reserve(krylov_dim);
        if (doPrecond_)
            precondVecSpace.reserve(krylov_dim);
        // vec_type y = vec_type::Zero(krylov_dim+1);
        vec_type y = vec_type(krylov_dim+1, 0.);
        std::vector<GivensType> givensRotations;
        givensRotations.reserve(krylov_dim);
        // mat_type H = mat_type::Zero(krylov_dim+1, krylov_dim);
        // mat_type R = mat_type::Zero(krylov_dim+1, krylov_dim);
        mat_type H(krylov_dim+1, krylov_dim, 0.);
        mat_type R(krylov_dim+1, krylov_dim, 0.);
        mps_tensor initialError = rhs - applyOperator(x), preconditionedError;
        double initialErrorNorm;
        if (doPrecond_) {
            preconditionedError = initialError;
            precond(preconditionedError);
            initialErrorNorm = std::sqrt(maquis::real(ietl::dot(initialError, preconditionedError)));
        }
        else {
            initialErrorNorm = ietl::two_norm(initialError);
        }
        double rhsNorm = ietl::two_norm(rhs);
        residual.push_back(initialErrorNorm);
        initialError /= initialErrorNorm;
        if (residual[0] > zeroThresh_) {
            if (doPrecond_) {
                preconditionedError /= initialErrorNorm;
                vecSpace.push_back(preconditionedError);
                precondVecSpace.push_back(initialError);
            }
            else {
                vecSpace.push_back(initialError); 
            }
        }
        y[0] = residual[0];
        while (residual[iter] > tol && iter < krylov_dim-1 && !exit) {
            maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                         << residual[iter] << std::endl;
            // Begin of the Arnoldi part
            auto Av = applyOperator(vecSpace[iter]);
            if (iter > 0) {
                H(iter-1, iter) = H(iter, iter-1);
                if (doPrecond_)
                    Av -= H(iter, iter-1)*precondVecSpace[iter-1];
                else
                    Av -= H(iter, iter-1)*vecSpace[iter-1];
            }
            auto alpha = ietl::dot(vecSpace[iter], Av);
            if (doPrecond_)
                Av -= alpha*precondVecSpace[iter];
            else
                Av -= alpha*vecSpace[iter];
            H(iter, iter) += alpha;
            mps_tensor pAv;
            if (doPrecond_) {
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
                if (doPrecond_) {
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
            residual.push_back(std::abs(y[iter+1])/rhsNorm);
            iter += 1;
        }
        // Final back-substitution
        if (iter != 0) {
            //vec_type result = R.block(0, 0, iter, iter).colPivHouseholderQr().solve(y.head(iter));
            //vec_type result = vec_type::Zero(iter);
            vec_type result(iter, 0.);
            for (int i = 0; i < iter; i++) {
                scalar_type s = 0.;
                for (int j = i+1; j < iter; j++)
                    s = s + R(i, j)*result[j];
                result[i] = (y[i] - s)/R(i, i);
            }
            for (int iFinal = 0; iFinal < iter; iFinal++)
                x += result[iFinal]*vecSpace[iFinal];
        }
        printEndl();
    }
    
    /**
     * @brief Solve the local linear system with the MINRES algorithm.
     *
     * The implementation is based on the PyKry python library. 
     *
     * @param sp site problem (representing the Hamiltotnian.
     * @param rhs right-hand side of the linear system.
     * @param krylov_dim maximum size of the Krylov subspace.
     * @param tol convergence tolerance for the whole GMRES procedure
     */
    void minres(const SiteProblem<Matrix, SymmGroup>& sp, const mps_tensor& rhs,
                int krylov_dim, real_type tol) {
        // Printing
        printHeader();
        maquis::cout << " --------------------------------- " << std::endl;
        maquis::cout << " Iteration  | Rel. error estimate  " << std::endl;
        maquis::cout << " --------------------------------- " << std::endl;
        // Sets up the initial value of all parameters.
        int iter = 0;
        bool exit = false;
        std::vector<mps_tensor> vecSpace, precondVecSpace;
        std::vector<double> residual;
        residual.reserve(krylov_dim);
        vecSpace.reserve(krylov_dim);
        if (doPrecond_)
            precondVecSpace.reserve(krylov_dim);
        // mat_type H = mat_type::Zero(krylov_dim+1, krylov_dim);
        mat_type H(krylov_dim+1, krylov_dim, 0.);
        auto p = rhs;
        auto rhsNorm = ietl::two_norm(rhs);
        auto v = applyOperator(p);
        residual.push_back(ietl::two_norm(rhs - v));
        scalar_type vnorm = ietl::dot(p, v);
        if (std::norm(vnorm) > zeroThresh_) {
            precondVecSpace.push_back(p/vnorm);
            vecSpace.push_back(v/vnorm);
        }
        auto W = std::vector<mps_tensor>(2);
        std::vector<scalar_type> y(2);
        y[0] = residual[0];
        auto yk = 0.*rhs;
        GivensType G1, G2;
        mps_tensor MAv;
        while (residual[iter] > tol && iter < krylov_dim-1 && !exit) {
            maquis::cout << std::setw(5) << iter << "          " << std::setw(15) << std::scientific
                         << residual[iter] << std::endl;
            auto Av = applyOperator(vecSpace[iter]);
            if (iter > 0) {
                H(iter-1, iter) = H(iter, iter-1);
                if (doPrecond_) 
                    Av -= H(iter, iter-1)*precondVecSpace[iter-1];
                else
                    Av -= H(iter, iter-1)*vecSpace[iter-1];
            }
            auto alpha = ietl::dot(vecSpace[iter], Av);
            if (doPrecond_) 
                Av -= alpha*precondVecSpace[iter];
            else
                Av -= alpha*vecSpace[iter];
            H(iter, iter) = alpha;
            if (doPrecond_) {
                MAv = Av;
                precond(MAv);
                H(iter+1, iter) = ietl::dot(Av, MAv);
            } else {
                H(iter+1, iter) = ietl::dot(Av, Av);
            }
            // Check if it's zero norm
            if (std::abs(H(iter+1, iter)) < zeroThresh_) {
                exit = true; 
            } 
            else {
                if (doPrecond_) {
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
            R(2) = std::real(H(iter, iter));
            R(3) = std::real(H(iter+1, iter));
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
            residual.push_back(std::abs(y[0])/rhsNorm);
        }
        x += yk;
        printEndl();
    }

public:
    /**
     * @brief Apply the shifted operator onto the MPS.
     * @param inputVec lhs vector
     * @return mps_tensor rhs vector
     */
    mps_tensor applyOperator(const mps_tensor& inputVec) const {
        mps_tensor ret;//, retSquared;
        ietl::mult(sp, inputVec, ret);
        //if (isFolded_) {
        //    ietl::mult(sp, inputVec, retSquared, 0, true);
        //    if (params["lin_alg"] == "feast")
        //        ret = retSquared - 2.*maquis::real(Z)*ret + inputVec*std::norm(Z);
        //    else
        //        ret = retSquared - 2.*sigma*ret + inputVec*sigma*sigma;
        //}
        //else {
        if (params["linsystem_dmrg_alg"] == "feast")
            ret = ret - Z*inputVec;
        else
            ret = ret - sigma*inputVec;
        //}
        return ret;
    }

    /** @brief Overloading of user-defined conjugate function */
    static inline double localConj(double a) {
        return a;
    }
    
    /** @brief Conjugate function for complex numbers */
    static inline std::complex<double> localConj(std::complex<double> a) {
        return std::conj(a);
    }

private:

    /** @brief Preconditioner */
    void precond(mps_tensor& inputVec) const {
        block_matrix<Matrix, SymmGroup>& data = inputVec.data();
        scalar_type denom;
        assert(shape_equal(data, precond_));
        for (size_t b = 0; b < data.n_blocks(); ++b) {
            for (size_t i = 0; i < num_rows(data[b]); ++i) {
                for (size_t j = 0; j < num_cols(data[b]); ++j) {
                    if (params["lin_alg"] == "feast")
                        denom = (precond_[b](i, j) - Z);
                    else
                        denom = (precond_[b](i, j) - sigma);
                    if (std::fabs(denom) > 1.0E-10)
                        data[b](i, j) /= std::fabs(denom);
                }
            }
        }
    }

    /** @brief (complex-valued) Givens rotation (needed for MINRES) */ 
    std::tuple<double, scalar_type, scalar_type> symGivens(scalar_type a, scalar_type b) {
        double c;
        scalar_type s, r;
        if (std::abs(b) < 1.0E-16) {
            c = 1;
            s = static_cast<scalar_type>(0.);
            r = a;
        }
        else if (std::abs(a) < 1.0E-16) {
            c = 0;
            s = static_cast<scalar_type>(1.);
            r = b;
        }
        else if (std::abs(b) > std::abs(a)) {
            auto t = std::abs(a) / std::abs(b);
            c = 1./std::sqrt(1. + std::norm(t));
            s = c * localConj((b/std::abs(b)) / (a/std::abs(a)));
            c = c*t;
            r = b / localConj(s);
        }
        else {
            auto t = std::abs(b) / std::abs(a);
            c = 1./std::sqrt(1. + std::norm(t));
            s = c * t * localConj( (b/std::abs(b)) / (a/std::abs(a)));
            r = a / c;
        }
        return std::make_tuple(c, s, r);
    }

    /** @brief Just prints a line for the table of the results */
    void printEndl() const {
        maquis::cout << " --------------------------------- " << std::endl;
        maquis::cout << std::endl;
        maquis::cout << std::fixed;
    }

    /** @brief Prints the header of the table */
    void printHeader() const {
        maquis::cout << std::endl;
        maquis::cout << " +-----------------------------------------+ " << std::endl;
        maquis::cout << " | ITERATIVE SOLUTION OF THE LINEAR SYSTEM | " << std::endl;
        maquis::cout << " +-----------------------------------------+ " << std::endl;
    }

    /** @brief Calculates the Givens rotation */
    std::tuple<scalar_type, scalar_type, scalar_type> givens(scalar_type x0, scalar_type x1, scalar_type in0, scalar_type in1) {
        double cos;
        scalar_type sin;
        std::tie(cos, sin) = calculateRotation(x0, x1);
        auto r = cos*x0 + sin*x1;
        auto out0 = cos*in0 + sin*in1;
        auto out1 = -localConj(sin)*in0+cos*in1;
        return std::make_tuple(out0, out1, r);
    }


    /* Private members */
    SiteProblem<Matrix, SymmGroup> const& sp;
    block_matrix<Matrix, SymmGroup> const& precond_;
    //mps_tensor const & rhs;
    mps_tensor rhs;
    BaseParameters& params;
    double sigma;
    scalar_type Z;
    energy_type en;
    mps_tensor x;
    bool doPrecond_, isFolded_;
    static constexpr double zeroThresh_ = 1.0E-16;
};

#endif
