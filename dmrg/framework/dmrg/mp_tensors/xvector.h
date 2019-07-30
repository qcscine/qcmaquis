/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#ifndef XVECTOR_H
#define XVECTOR_H

#include "mpstensor.h"
#include "mps.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

namespace lr {

    // Functions required in the calculation of nonredundant MPS parameters for the linear response

    // Class to store nonredundant MPS parameters.
    // Must be a class rather than a simple std::vector of block matrices
    // because we must store the first parameter (M) as an MPSTensor with all its internal properties
    // and also some transformation matrices
    template <class Matrix, class SymmGroup>
    class XVector
    {
        public:
        // choose abelian or SU2 GEMM depending on SymmGroup
        // typedef typename boost::mpl::if_<symm_traits::HasSU2<SymmGroup>, SU2::SU2Gemms, contraction::abelian::Gemms>::type::gemm gemm;

        typedef typename alps::numeric::associated_real_diagonal_matrix<Matrix>::type DiagMatrix;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > bm_vector;

        // Iterators for range-based access
        typedef typename bm_vector::iterator iterator;
        typedef typename bm_vector::const_iterator const_iterator;

        iterator begin() { return X_.begin(); };
        const_iterator begin() const { return X_.begin(); };
        iterator end() { return X_.end(); };
        const_iterator end() const { return X_.end(); };

        // Empty constructors
        XVector() : length_(0), X_(), S_(), N_() {};
        XVector(std::size_t length) : length_(length), X_(length), S_(length-1), N_(length) {};

        // Constructor with an implicit conversion from B to X

        // From a complete right-canonized MPS, obtain a vector with nonredundant MPS parameters, the
        // MPSTensor for the first site and the X matrices for all remaining sites with
        //                  3m      m?                                  sigma_i          sigma_i
        // X             = sum     sum    sum      S^(1/2)             B               N+
        //  a_i-1 a'_i-1   a''_i-1  a''_i sigma_i        a_i-1 a''_i-1  a''_i-1 a''_i    a''_i a'_i-1

        // We need to work with a copy of the MPS, since we will shift canonisation
        XVector(MPS<Matrix, SymmGroup> mps) : XVector(mps.length())
        {
            // just to make sure the MPS is right-normalised
            mps.normalize_right();

            assert(length_ > 0);

            // construct X from B
            for (std::size_t i = 0; i < length_; i++)
                {
                    block_matrix<Matrix, SymmGroup> U; // eigenvectors
                    block_matrix<Matrix, SymmGroup> tmp;

                    // Build the overlap matrix only for i > 0, for site 0 it should be unity
                    if (i > 0)
                    {
                        // Build the overlap matrix for the right renormalised states
                        // i.e. build MM^T from right-paired MPS tensors at site i-1
                        // dm            = sum              M                    M^T
                        //   a_i-1 a_i-1'  (a_i sigma_i)     a_i-1 (a_i sigma_i)  (a_i sigma_i) a_i-1'

                        mps.canonize(i);
                        mps[i].make_right_paired();
                        block_matrix<Matrix, SymmGroup> dm;
                        gemm(mps[i].data(), transpose(conjugate(mps[i].data())), dm);

                        // build S=dm^(1/2)
                        block_matrix<DiagMatrix, SymmGroup> Lambda; // eigenvalues
                        heev(dm, U, Lambda);
                        // zero_small_values_inplace(Lambda); // set very small negative values to zero because otherwise sqrt does not work // should never be needed because we must be able to invert Lambda
                        const block_matrix<DiagMatrix, SymmGroup>& sqrtLambda = sqrt(Lambda);
                        gemm(U,sqrtLambda,tmp);
                        gemm(tmp,transpose(conjugate(U)), S_[i-1]);
                    }

                    // Build N
                    N_[i] = constructN(mps[i]);

                    // Construct N^+ from N with a pseudoinverse using SVD
                    block_matrix<Matrix, SymmGroup> V;
                    block_matrix<DiagMatrix, SymmGroup> sigma;
                    // SVD of N: N=USV
                    svd(N_[i], U, V, sigma);

                    block_matrix<Matrix, SymmGroup> Ninv;
                    // Build the pseudoinverse: N+=V^T S^(-1) U^T
                    // block_matrix<Matrix, SymmGroup> invSigma = inverse(sigma); // Workaround because SU2 gemm doesn't work with diagonal matrices
                    const block_matrix<DiagMatrix, SymmGroup> & invSigma = inverse(sigma);
                    block_matrix<Matrix, SymmGroup> VT = transpose(conjugate(V));
                    gemm(VT,invSigma,tmp);
                    gemm(tmp,transpose(conjugate(U)), Ninv);

                    // Calculate the pseudoinverse as N^T *(NN^T)^-1 to check if the pseudoinverse with SVD yields the correct result

                    // block_matrix<Matrix, SymmGroup> NNT, NNT_U, tmpNNT,tmpNNT2,Ninv2;
                    // block_matrix<DiagMatrix, SymmGroup> NNT_Lambda;
                    // // NN^T
                    // gemm(N_[i-1],transpose(conjugate(N_[i-1])),NNT);
                    // heev(NNT, NNT_U, NNT_Lambda);
                    // //(NN^T)^-1
                    // gemm(NNT_U,inverse(NNT_Lambda),tmpNNT);
                    // gemm(tmpNNT,transpose(conjugate(NNT_U)),tmpNNT2);
                    // //N^T(NN^T)^-1
                    // gemm(transpose(conjugate(N_[i-1])), tmpNNT2, Ninv2);
                    // It seems so, since Ninv = Ninv2

                    // Now that we have N^+ (in Ninv) and S^(1/2) (in S), contract the matrices to obtain X
                    // Multiply B with N^+: B must be right-paired but Ninv (probably?) left-paired
                    mps[i].make_right_paired();

                    // Premultiply the result with S and store it in X
                    // For site 0 S is unity so ignore it
                    if (i == 0)
                        gemm(mps[i].data(), Ninv, this->X_[i]);
                    else
                    {
                        gemm(mps[i].data(), Ninv, tmp);
                        gemm(S_[i-1], tmp, this->X_[i]);
                    }


                }
        };

        // Move constructor
        XVector(XVector&& other) : XVector() // initialize via default constructor
        {
            swap(*this, other);
        }

        // Assignment operator
        XVector& operator=(XVector other)
        {
            swap(*this, other);
            return *this;
        }

        // Get a right-normalised MPSTensor B formed from nonredundant MPS parameters X at site i
        // Transformation from X to B
        // From a set of non-redundant MPS consisting of an MPSTensor M for the first site and X matrices
        // obtain a right-canonized MPS, i.e. transform all X to right canonised MPS tensors B with
        //  sigma_i      m      3m                                              sigma_i
        // B          = sum     sum     S^(-1/2)               X               N
        //  a_i-1 a_i   a'_i-1  a''_i-1         a_i-1 a'_i-1    a'_i-1 a''_i-1  a''_i-1 a_i

        MPSTensor<Matrix,SymmGroup> getB(std::size_t site) const
        {
            MPSTensor<Matrix, SymmGroup> ret;
            block_matrix<Matrix, SymmGroup> B;

            if (site == 0)
            {
                // perform the transformation from X to B above
                // For site 0, do not multiply with S as it is unity
                gemm(X_[0], N_[0], B);
                ret.replace_right_paired(B);
            }
            else
            {
                //Build S^(-1/2) from S^1/2
                block_matrix<Matrix, SymmGroup> U, tmp, S;
                block_matrix<DiagMatrix, SymmGroup> Lambda;
                heev(S_[site-1], U, Lambda);

                const block_matrix<DiagMatrix, SymmGroup>& invLambda = inverse(Lambda);
                gemm(U, invLambda, tmp);
                // S contains S^(-1/2) now
                gemm(tmp,transpose(conjugate(U)), S);

                // Contract S with X and N to get B as block_matrix
                gemm(X_[site], N_[site], tmp);
                gemm(S,tmp, B);

                // form a right-paired MPS tensor from B
                ret.replace_right_paired(B);
            }
            return ret;
        }

        MPS<Matrix, SymmGroup> transformXtoB() const
        {
            assert(length_ > 0);
            MPS<Matrix,SymmGroup> ret(length_);

            for (std::size_t i = 0; i < length_; i++)
                ret[i] = getB(i);
            return ret;
        }

        // return either an MPSTensor for site 0 or block_matrix for other sites
        // also include an out of bounds check (slow?)
        block_matrix<Matrix, SymmGroup>& operator[](std::size_t i)
        {
            return X_[i];
        }

        // loading and saving of XVectors from files
        // shamelessly copy-pasted from mps.hpp

        void load(std::string const& dirname)
        {
            /// get size of MPS
            std::size_t L = 0;
            while (boost::filesystem::exists( dirname + "/xvec" + boost::lexical_cast<std::string>(++L) + ".h5" ));

            /// load tensors
            XVector<Matrix, SymmGroup> tmp(L);
            size_t loop_max = tmp.length();
            parallel::scheduler_balanced scheduler(loop_max);
            for(size_t k = 0; k < loop_max; ++k){
                parallel::guard proc(scheduler(k));
                std::string fname = dirname+"/xvec"+boost::lexical_cast<std::string>((size_t)k)+".h5";
                storage::archive ar(fname);
                ar["/tensor"] >> tmp[k];
            }
            swap(*this, tmp);
        }

        void save(std::string const& dirname)
        {
            /// create chkp dir
            if(parallel::master() && !boost::filesystem::exists(dirname))
                boost::filesystem::create_directory(dirname);

            parallel::scheduler_balanced scheduler(length_);
            size_t loop_max = length_;

            for(size_t k = 0; k < loop_max; ++k){
                parallel::guard proc(scheduler(k));
                storage::migrate(X_[k]);
            }
            parallel::sync();

            for(size_t k = 0; k < loop_max; ++k){
                parallel::guard proc(scheduler(k));
                if(!parallel::local()) continue;
                const std::string fname = dirname+"/xvec"+boost::lexical_cast<std::string>((size_t)k)+".h5.new";
                storage::archive ar(fname, "w");
                ar["/tensor"] << X_[k];
            }

            parallel::sync(); // be sure that chkp is in valid state before overwriting the old one.

            omp_for(size_t k, parallel::range<size_t>(0,loop_max), {
                parallel::guard proc(scheduler(k));
                if(!parallel::local()) continue;
                const std::string fname = dirname+"/xvec"+boost::lexical_cast<std::string>((size_t)k)+".h5";
                boost::filesystem::rename(fname+".new", fname);
            });
        }

        std::size_t length() const { return length_; };

        friend void swap(XVector& a, XVector& b)
        {
            using std::swap;
            swap(a.X_, b.X_);
            swap(a.S_, b.S_);
            swap(a.N_, b.N_);
            swap(a.length_, b.length_);
        }

        private:
            // Number of sites
            std::size_t length_;

            // MPS parameters
            bm_vector X_;

            // Transformation matrices N and S^-1/2 that must be stored for the back transformation
            bm_vector S_;
            bm_vector N_;

            // For a given MPS constructs the N matrix, which is the transformation matrix of the |a_i-1 sigma_i(R)>
            // basis to the basis COMPLEMENTARY to |a_i(R)> basis. (In contrast, in a regular DMRG sweep we
            // construct a transformation matrix from the |a_i-1 sigma_i(R)> to the |a_i(R)> basis). Its size is
            // therefore 3m x 4m.
            // It is assumed that all MPSTensors right of M are right-normalised.
            block_matrix<Matrix, SymmGroup> constructN(const MPSTensor<Matrix, SymmGroup> & M)
            {

                // Build the density matrix for the right renormalised states

                M.make_right_paired();
                block_matrix<Matrix, SymmGroup> dm;
                gemm(transpose(conjugate(M.data())), M.data(), dm);
                // Diagonalise the RDM.

                block_matrix<Matrix, SymmGroup> U; // eigenvectors
                block_matrix<DiagMatrix, SymmGroup> Lambda; // eigenvalues
                heev(dm, U, Lambda);

                // Remove the eigenvalues and the eigenvectors corresponding to the |a_i(R)> basis.

                const std::vector<std::size_t> & left_basis_sizes = M.data().left_basis().sizes();
                assert(left_basis_sizes.size() == U.n_blocks());
                for (std::size_t i = U.n_blocks()-1; i >= 0; i--)
                {
                    // If the size of the right basis is equal to that of the left basis (which shouldn't happen, right?)
                    if (U.right_basis().sizes()[i] <= left_basis_sizes[i])
                        U.remove_block(U.basis().left_charge(i),
                                    U.basis().right_charge(i));
                    else
                        // remove columns that are eigenvectors corresponding to the |a_i> basis.
                        U.remove_cols_from_block(i, 0, left_basis_sizes[i]);
                }

                U.transpose_inplace();
                // N [U] has 3 indices, just as an MPS tensor. If we were to save it in an MPSTensor, it would be a right-paired tensor
                return U;
            }
    };

    // Transformation from B to X implemented in the XVector class constructor
    template <class Matrix, class SymmGroup>
    XVector<Matrix, SymmGroup> transformBtoX(const MPS<Matrix, SymmGroup> & mps)
    {
        return XVector<Matrix, SymmGroup>(mps);
    }

};
#endif