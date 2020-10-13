// jdqr() according to http://web.eecs.utk.edu/~dongarra/etemplates/node144.html
/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2011 by Robin Jaeger <jaegerr@phys.ethz.ch>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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
#ifndef JACOBI_DAVIDSON_H
#define JACOBI_DAVIDSON_H
#include <alps/config.h> // needed to set up correct bindings

#include <ietl/traits.h>
#include <ietl/complex.h>
#include <ietl/iteration.h>
#include <ietl/vectorspace.h>
#include <ietl/gmres.h>
#include <ietl/bicgstabl.h>
#include <vector>
#include <cassert>
#include <algorithm>
#include <exception>
#include <limits>

#include <boost/function.hpp>

#include <boost/numeric/ublas/hermitian.hpp>
#include <boost/numeric/ublas/matrix.hpp>
#include <boost/numeric/ublas/matrix_proxy.hpp>
#include <boost/numeric/ublas/io.hpp>
//for lapack::heev
#include <boost/numeric/bindings/lapack/driver/heev.hpp>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
//for lapack::getrs
#include <boost/numeric/bindings/lapack/computational.hpp>
//#include <boost/numeric/bindings/detail/config/fortran.hpp>
#include <boost/numeric/bindings/ublas/matrix.hpp>
#include <boost/numeric/bindings/ublas/vector.hpp>
#include <boost/numeric/bindings/ublas/hermitian.hpp>
#include <boost/numeric/bindings/std.hpp>
#include <boost/numeric/bindings/ublas.hpp>
#include <boost/numeric/bindings/upper.hpp>
#include <boost/numeric/bindings/lower.hpp>


namespace ietl{

    namespace ublas=boost::numeric::ublas;
    namespace lapack=boost::numeric::bindings::lapack;

    
//jd_iterator//////////////////////////////////////////////////////
    template <class T>
    class jd_iteration : public basic_iteration<T> {
        public:
        jd_iteration(   size_t max_iter, 
                        size_t m_min = 10, 
                        size_t m_max = 20, 
                        T reltol = sqrt(std::numeric_limits<double>::epsilon()), 
                        T abstol = sqrt(std::numeric_limits<double>::epsilon())
                    ) 
            : basic_iteration<T>(max_iter, reltol, abstol), m_min_(m_min), m_max_(m_max) {}

        inline size_t m_min() const
            {   return m_min_;  }
        inline size_t m_max() const
            {   return m_max_;  }
        const std::string& error_msg() 
            { return basic_iteration<T>::err_msg; }    
        private:
        size_t m_min_, m_max_;
    };
//jd_iterator//////////////////////////////////////////////////////

    namespace detail{

//deflated matrix multiplication (I-QQ*)(A-\theta I)(I-QQ*)////////

        template <class MATRIX, class VS> 
            class deflated_matrix;
        template <class MATRIX, class VS, class VECTOR>
            void mult(const deflated_matrix<MATRIX,VS> & Adef, const VECTOR& v, VECTOR& r);

        template <class MATRIX, class VS>
        class deflated_matrix   {
            public:
                typedef typename vectorspace_traits<VS>::vector_type vector_type;
                typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
                typedef typename real_type<scalar_type>::type real_type;

                deflated_matrix(const MATRIX& A, double theta, const std::vector<vector_type>& Q) 
                    : A_(A) , theta_(theta) , Q_(Q) {}
                ~deflated_matrix() {}

                void set_theta(real_type theta) 
                    { theta_ = theta; }

                friend void mult<>(const deflated_matrix& A_def, const vector_type& v, vector_type& r);

            private:
                const MATRIX& A_;
                real_type theta_;
                const std::vector<vector_type>& Q_;    
        };

        template <class MATRIX, class VS, class VECTOR>
        void mult(  const deflated_matrix<MATRIX,VS>& Adef, 
                    const VECTOR& v, 
                    VECTOR& r)
        {// (I-QQ*)(A-thetaI)(I-QQ*)v = r
                VECTOR temp = v;

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    temp -= ietl::dot(Adef.Q_[i], temp) * Adef.Q_[i];

                ietl::mult(Adef.A_, temp, r);
                r -= Adef.theta_ * temp;

                for(size_t i = 0; i < Adef.Q_.size() ; ++i) 
                    r -= ietl::dot(Adef.Q_[i], r) * Adef.Q_[i];
        }
        
        //vector-set - matrix multiplication, you may want to provide your own
        //for your specific vector type
        //overwrites vecset
        template <class VECTOR, class MATRIX>
        void mult( std::vector<VECTOR>& vecset, const MATRIX& mat)
        {   // mat.size1() = m, mat.size2() = m-1
            assert(vecset.size() == mat.size1());
            std::vector<VECTOR> tmp(mat.size2());
            for(std::size_t i = 0; i < tmp.size(); ++i)
            {
                tmp[i] = vecset[0] * mat(0,i);
                for(std::size_t j = 1; j < vecset.size(); ++j)
                    tmp[i] += vecset[j] * mat(j,i);
            }
            std::swap(vecset, tmp);
        }
        //vector-set - matrix multiplication for boost::numeric::ublas::vector<T>
        template <class T, class MATRIX>
        void mult( std::vector< ublas::vector<T> >& vecset, const MATRIX& mat)
        {
            assert(vecset.size() == mat.size1());
            for (std::size_t i = 0; i < vecset[0].size(); ++i) 
            {
                ublas::vector<T> tmp(mat.size2(), 0);

                for (std::size_t j = 0; j < mat.size2(); ++j)
                    for (std::size_t k = 0; k < vecset.size(); ++k)
                        tmp(j) += vecset[k](i) * mat(k, j);

                for (std::size_t j = 0; j < mat.size2(); ++j)
                    vecset[j](i) = tmp(j);
            }

            vecset.resize(mat.size2());
        }

}//end namespace detail///////////////////////////////////////////////////////
//solver//////////////////////////////////////////////////////////////////////
namespace solver {

//left preconditioned correction-equation-solver
    template <class SOLV, class MATRIX, class VS, class PREC>
        class left_prec_solver;

    template <class SOLV, class MATRIX, class VS, class PREC, class VECTOR>
    void mult (const left_prec_solver<SOLV,MATRIX,VS,PREC>&, const VECTOR&, VECTOR&);

    template <class SOLV, class MATRIX, class VS, class PREC>
    class left_prec_solver {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        left_prec_solver (SOLV& s, const MATRIX& A, const std::vector<vector_type>& Q, const VS& vs, PREC& K) 
            : solv_(s), A_(A) , Q_(Q) , vspace_(vs) , K_(K) , x0_(new_vector(vs))  {}

        template <class IT>
        void operator() ( real_type theta, vector_type& r, vector_type& t, IT& iter );

        friend void mult<>(const left_prec_solver& A_def, const vector_type& v, vector_type& z);

    protected:
        SOLV& solv_;
        const MATRIX& A_;
        const std::vector<vector_type>& Q_;
        const VS& vspace_;
        std::vector<vector_type> Q_hat_;
        PREC& K_;
        vector_type x0_;
        real_type theta_;
        mutable ublas::vector<scalar_type> gamma_;
        std::vector<fortran_int_t> pivot_;
        matrix_type LU_, M_;
    };

    template <class SOLV, class MATRIX, class VS, class PREC, class VECTOR>
    void mult (const left_prec_solver<SOLV,MATRIX,VS,PREC>& LPS, const VECTOR& v, VECTOR& z)
    {
        VECTOR temp;

        z = LPS.theta_ * v;

        ietl::mult(LPS.A_, v, temp);

        temp -= z;
        // K y_tilde = y
        ietl::mult(LPS.K_, temp, z);

        for(size_t i = 0; i < LPS.Q_.size(); ++i)
            LPS.gamma_(i) = ietl::dot(LPS.Q_[i],z);

        int info = lapack::getrs(LPS.LU_, LPS.pivot_, LPS.gamma_);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        for(size_t i = 0; i < LPS.Q_hat_.size(); ++i)
            z -= LPS.Q_hat_[i]* LPS.gamma_(i);
    }

    // provide a function ietl::mult(Prec K, v1, v2) with v2 ~= K v1
    template <class SOLV, class MATRIX, class VS, class PREC> template <class IT>
    void left_prec_solver<SOLV,MATRIX,VS,PREC>::operator() ( real_type theta, vector_type& r, vector_type& t, IT& iter )
    {
        theta_ = theta;
        const fortran_int_t m = Q_.size(), m_old = Q_hat_.size();//if K changes =0;
        double abs_tol = iter.absolute_tolerance();
        abs_tol *= abs_tol;

        M_.resize(m,m);
        gamma_.resize(m);
        pivot_.resize(m);
        Q_hat_.resize(m);

        for(size_t i = m_old; i < m; ++i)
            ietl::mult(K_, Q_[i], Q_hat_[i]);

        for(size_t i = 0; i < m; ++i)
            for(size_t j = ( (i<m_old) ? m_old : 0 ); j < m; ++j)
                M_(i,j) = ietl::dot(Q_[i], Q_hat_[j]);

        LU_ = M_;

        int info = lapack::getrf(LU_, pivot_);
            if(info != 0)   throw std::runtime_error("lapack::getrf failed.");

        ietl::mult(K_, r, t);

        for(size_t i = 0; i < m; ++i)
            gamma_(i) = ietl::dot(Q_[i],t);

        //calculate alpha from M alpha = gamma
        info = lapack::getrs(LU_, pivot_, gamma_);
            if(info != 0)   throw std::runtime_error("lapack::getrs failed.");

        r = -t;

        for(size_t i = 0; i < m; ++i)
            r += gamma_(i)*Q_hat_[i];

        t = solv_ (*this, r, x0_, abs_tol);

    }//left_prec_solver::void()


// non-preconditioned simple solver//////////////////////////////////////////////////////
    template <class SOLV, class MATRIX, class VS>
    class jd_solver {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename real_type<scalar_type>::type real_type;

        jd_solver (SOLV& s, const MATRIX& A, const std::vector<vector_type>& Q, const VS& vspace)
            : Adef_(A, 0, Q), vspace_(vspace), x0_(new_vector(vspace)), Q_(Q), solv_(s)
             {}
        
        template <class IT>
        void operator() ( real_type theta, vector_type& r, vector_type& t, IT& iter );

    protected:
        detail::deflated_matrix<MATRIX, VS> Adef_;
        const VS& vspace_;
        vector_type x0_;
        const std::vector<vector_type>& Q_;
        SOLV& solv_;
    };
        template <class SOLV, class MATRIX, class VS>
        template <class IT>
        void jd_solver<SOLV,MATRIX,VS>::operator()( real_type theta, vector_type& r, vector_type& t, IT& iter)
        {
            Adef_.set_theta(theta);

            //one step approximation
            /*x0_ = -r;
            for(size_t i = 0; i < Q_.size(); ++i)
                x0_ += ietl::dot(Q_[i],r)/ietl::dot(Q_[i],Q_[i]) * Q_[i];*/

            r *= -1;

            //starting with x0_ = 0
            t = solv_ (Adef_, r, x0_, iter.absolute_tolerance()*iter.absolute_tolerance());
        }
}//end namespace solver///////////////////////////////////////////////////////////

    template <class MATRIX, class VS>
    class jd {
    public:
        typedef typename vectorspace_traits<VS>::vector_type vector_type;
        typedef typename vectorspace_traits<VS>::scalar_type scalar_type;
        typedef typename vectorspace_traits<VS>::size_type size_type;
        typedef typename vectorspace_traits<VS>::magnitude_type magnitude_type;
        typedef typename real_type<scalar_type>::type real_type;
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;

        jd (const MATRIX& A, VS& vspace, size_t v = 0)
            : A_(A), vspace_(vspace), n_(vec_dimension(vspace)), verbose_(v) {}

        //default: search for lowest eigenvalues


        //without preconditioning
        #ifndef __GXX_EXPERIMENTAL_CXX0X__
        template <class SOLV, class IT, class GEN>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, SOLV f, bool search_highest = false)
        #else
        template <class SOLV=ietl_gmres, class IT, class GEN>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, SOLV f = ietl_gmres(), bool search_highest = false)
        #endif
            {
                typedef solver::jd_solver<SOLV,MATRIX,VS> solver_type;
                solver_type solver(f, A_, X_, vspace_);
                jdqr <solver_type,IT,GEN> (solver, iter, gen, k_max, search_highest);
            }

        //with left preconditioner K with K ~= inv (A - \theta I)
        template <class SOLV, class IT, class GEN, class PREC>
        void eigensystem(IT& iter, GEN& gen, size_type k_max, PREC& K, SOLV f, bool search_highest = false)
            {
                typedef solver::left_prec_solver<SOLV,MATRIX,VS,PREC> solver_type;
                solver_type lp_solver(f, A_, X_, vspace_, K);
                jdqr <solver_type,IT,GEN> (lp_solver, iter, gen, k_max, search_highest);
            }

        //jdqr with harmonic ritz values
        template <class SOLV, class IT, class GEN>
        void eigensystem_harmonic(IT& iter, GEN& gen, size_type k_max, real_type tau, SOLV& f)    
            {
                typedef solver::jd_solver<SOLV,MATRIX,VS> solver_type;
                solver_type solver(f, A_, X_, vspace_);
                jdqr_harmonic <solver_type,IT,GEN> (solver, iter, gen, k_max, tau);
            }

        //access functions
        std::pair<real_type,vector_type> eigenpair(size_type k)
            {
                assert (k < Lambda_.size());
                return std::make_pair (Lambda_[k], X_[k]);
            }
        real_type eigenvalue(size_type k)
            {
                assert (k < Lambda_.size());
                return Lambda_[k];
            }
        real_set_type eigenvalues()
            {
                return Lambda_;
            }
        vector_type eigenvector (size_type k)
            {
                assert (k < X_.size());
                return X_[k];
            }
        void reset()
            {
                X_.clear();
                Lambda_.clear();
            }

    protected:
        const MATRIX& A_;
        VS& vspace_;
        size_type n_;
        vector_set_type X_;
        real_set_type Lambda_;
        size_t verbose_;    //verbosity levels: 0 = no, 1 = yes, 2 = all
        
        template <class SOLVER, class IT, class GEN>
        void jdqr(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, bool search_highest);
        
        template <class SOLVER, class IT, class GEN>
        void jdqr_harmonic(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau);
    };// end of class jd

//jdqr for exterior eigenvalues ///////////////////////////////////////////////////
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN>
    void jd<MATRIX,VS>::jdqr(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, bool search_highest)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::upper > herm_matrix_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        assert(k_max <= n_);

        const magnitude_type kappa = 0.25; // kappa < 1
        magnitude_type norm_t, norm_r;
        const size_type m_min = iter.m_min(), m_max = iter.m_max();
        size_type k = 0, m = 0;
        int info;

        vector_set_type V;
        vector_set_type VA;
        V.reserve(m_max);
        VA.reserve(m_max);
        X_.reserve(X_.size()+k_max);
        
        X_.resize( X_.size()+1 ); //X_.back() == u

        herm_matrix_type M;

        real_set_type Theta;
        Theta.reserve(m_max);
        matrix_type S;

        vector_type uA, r;
        vector_type t = new_vector(vspace_);

        generate(t, gen);
        project(t, vspace_);

        //if eigenvectors already exist
        for(size_t i = 0; i < X_.size()-1;++i)
            t -= ietl::dot(t,X_[i])*X_[i];

        while(true) 
        {
            //modified Gram-Schmidt
            norm_t = ietl::two_norm(t);

            for(size_type i = 0; i < m; ++i) 
                t -= ietl::dot(V[i],t) * V[i];
            if(ietl::two_norm(t) <= kappa * norm_t)
                for(size_type i = 0; i < m; ++i)
                    t -= ietl::dot(V[i],t) * V[i];

            norm_t = ietl::two_norm(t);

            if(norm_t == 0/*< iter.absolute_tolerance()*/) {
                    std::cerr << "Correction vector t is zero:\n";
                    if(iter.first())
                        throw std::runtime_error("generating a starting vector failed.");
                    else
                        throw std::runtime_error("try to solve the correction equation more accurately.");
                    }

            V.push_back(t/norm_t);
            ++m;
            VA.resize(m);

            ietl::mult(A_, V[m-1], VA[m-1]);
            assert(m == V.size() && m == VA.size());

            M.resize(m,m);

            for(size_type i = 0; i < m; ++i)
                M(i, m-1) = ietl::dot(V[i], VA[m-1]);

            //eigendecomposition M = S THETA S*
            Theta.resize(m);

            S = M; //copy M

            info = lapack::heev('V', boost::numeric::bindings::upper(S), Theta); //eigenvectors in columns of S
                if(info < 0) throw std::runtime_error("lapack::heev - illegal value.");
                if(info > 0) throw std::runtime_error("lapack::heev - failed to converge.");
            
            //resort S if searching for highest eigenvalues
            if(search_highest)
            {
                std::reverse( Theta.begin(), Theta.end() );
                for(size_type i = 0; i < size_type(m/2); ++i)
                    swap( ublas::column(S, i), ublas::column(S, m-1-i) ); 
            }

            // u = V s_1
            X_.back() = V[0] * S(0,0);
            for(size_type j = 1; j < m; ++j) 
                X_.back() += V[j] * S(j,0);

            // uA = VA s_1
            uA = VA[0] * S(0,0);
            for(size_type j = 1; j < m; ++j) 
                uA += VA[j] * S(j,0);

            r = uA - Theta[0] * X_.back();

            // accept eigenpairs
            norm_r = ietl::two_norm(r);

            while(iter.finished(norm_r,Theta[0])) 
            {
            
                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << Theta[0] << "\tResidual: " << norm_r << "\n";

                Lambda_.push_back( Theta[0] );

                if(++k == k_max) return;

                if(m < 1)
                    throw std::runtime_error("search space is depleted, try a different generator.");

                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                ublas::matrix_range< matrix_type > Sproxy (S, ublas::range(0,S.size1()), ublas::range(1,S.size2()) );

                // v_i = V s_i+1
                detail::mult(V, Sproxy);
                
                // vA_i = VA s_i+1
                detail::mult(VA, Sproxy);
                
                Theta.erase( Theta.begin() );
                S.resize(m,m);         
                for(size_type i = 0; i < m; ++i)
                {
                    M(i,i) = Theta[i];
                    ublas::column(S, i) = ublas::unit_vector<scalar_type> (S.size1(), i);
                }
                assert ( S.size1() == m && S.size2() == m );

                X_.resize(X_.size()+1);
                X_.back() = V[0];
                r = VA[0] - Theta[0] * X_.back();
                norm_r = ietl::two_norm(r);

            } // accept eigenpairs

            // restart
            if(m >= m_max) {
                if(verbose_ > 1) 
                    std::cout<<"restart...\n";
                M = ublas::zero_matrix<scalar_type>(m_min,m_min);

                ublas::matrix_range< matrix_type > Sproxy (S, ublas::range(0,S.size1()), ublas::range(1,m_min) );

                detail::mult(V, Sproxy);
                    
                detail::mult(VA, Sproxy);

                V.insert(V.begin(), X_.back());
                ietl::mult(A_, X_.back(), uA);
                VA.insert(VA.begin(), uA);

                for(size_type i = 0; i < m_min; ++i)
                    M(i,i) = Theta[i];

                Theta.resize(m_min);
                m = m_min;
            }// restart

            //assure t is orthogonal to Q
            do {
              norm_t = ietl::two_norm(t);
              for(size_t i = 0; i< X_.size();++i)
                t -= ietl::dot(t,X_[i])*X_[i];
            } while (ietl::two_norm(t) < kappa * norm_t);

            // correction equation
            solver(Theta[0], r, t, iter); //solver is allowed to change r

            //assure t is orthogonal to Q
            do {
              norm_t = ietl::two_norm(t);
              for(size_t i = 0; i< X_.size();++i)
                t -= ietl::dot(t,X_[i])*X_[i];
            } while (ietl::two_norm(t) < kappa * norm_t);

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());

            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\t resid = " << norm_r << "\n";
        }// main loop

    }// jdqr for exterior evs

//jdqr for interior eigenvalues with harmonic ritz values////////////////////////////////////////
    template <class MATRIX, class VS>
    template <class SOLVER, class IT, class GEN >
    void jd<MATRIX,VS>::jdqr_harmonic(SOLVER& solver, IT& iter, GEN& gen, size_type k_max, real_type tau)
    {
        typedef std::vector<vector_type> vector_set_type;
        typedef std::vector<real_type> real_set_type;
        typedef ublas::hermitian_matrix< scalar_type, ublas::lower > herm_matrix_type;
        typedef ublas::matrix< scalar_type, ublas::column_major > matrix_type;

        herm_matrix_type M;
        matrix_type S;

        real_set_type Theta;

        vector_set_type V, W;
        V.reserve(iter.m_max());
        W.reserve(iter.m_max());
        vector_type t = new_vector(vspace_);

        X_.reserve(X_.size()+k_max);
        X_.resize(X_.size()+1);

        generate(t, gen);
        project(t, vspace_);

        //if eigenvectors already exist
        for(size_t i = 0; i < X_.size()-1;++i)
            t -= ietl::dot(t,X_[i])*X_[i];

        size_t m = 0, k = 0;

        while(true) //k<k_max
        {
            W.resize(m+1);
            ietl::mult(A_, t, W.back() );
            W.back() -= tau * t;

            for(size_t i=0; i < m; ++i)
            {
                scalar_type gamma = ietl::dot(W[i], W.back());
                W.back() -= gamma * W[i];
                t -= gamma * V[i];
            }

            ++m;
            real_type norm_w = ietl::two_norm(W.back());
            if(norm_w == 0)
                throw std::runtime_error("New search vector is zero.");

            W.back() /= norm_w;
            V.push_back( t/norm_w );

            M.resize(m,m);
            for(size_t i=0; i < m; ++i)
                M(i,m-1) = ietl::dot(W[i],V.back());
            
            //eigendecomposition M = S THETA S*
            Theta.resize(m);
            S = M; //copy M

            int info = lapack::heev('V', boost::numeric::bindings::upper(S), Theta); //eigenvectors in columns of S
                if(info < 0) throw std::runtime_error("lapack::heev - illegal value.");
                if(info > 0) throw std::runtime_error("lapack::heev - failed to converge.");
            //ascending instead of descending

            vector_type u_tilde = V[0] * S(0,m-1);
            for(size_t i=1; i < m; ++i)
                u_tilde += V[i] * S(i,m-1);

            magnitude_type mu = ietl::two_norm(u_tilde);

            X_.back() = u_tilde/mu;
            real_type theta = Theta.back()/mu/mu;

            vector_type w_tilde = W[0] * S(0,m-1);
            for(size_t i=1; i < m; ++i)
                w_tilde += W[i] * S(i,m-1);

            vector_type r = w_tilde/mu - theta * X_.back();

            magnitude_type norm_r = ietl::two_norm(r);

            while(iter.finished(norm_r,Theta[0])) 
            {
                if(iter.error_code() == 1)
                    throw std::runtime_error(iter.error_msg());

                if(verbose_ > 0)
                    std::cout << "Accepted eigenpair #"<< k+1 << "\t" << theta + tau << "\t Residual: " << norm_r << "\n";
                Lambda_.push_back( theta + tau );

                if(++k == k_max) return;

                if(m < 1)
                    throw std::runtime_error("search space is depleted, try a different generator.");
                --m;
                M = ublas::zero_matrix<scalar_type>(m,m);

                vector_set_type tempvs (m);
                for(size_t i = 0; i < m; ++i)
                {
                    tempvs[i] = V[0] * S(0,S.size2()-2-i);
                    for(size_t j = 1; j < S.size1(); ++j)
                        tempvs[i] += V[j] * S(j,S.size2()-2-i);
                }

                V.pop_back();
                for(size_t i = 0; i < m; ++i)
                    std::swap(V[i],tempvs[i]);
                for(size_t i = 0; i < m; ++i)
                {
                    tempvs[i] = W[0] * S(0,S.size2()-2-i);
                    for(size_t j = 1; j < S.size1(); ++j)
                        tempvs[i] += W[j] * S(j,S.size2()-2-i);
                }

                W.pop_back();
                for(size_t i = 0; i < m; ++i)
                    std::swap(W[i],tempvs[i]);
                Theta.pop_back();
                S.resize(m,m);
                for(size_t i = 0; i < m; ++i)
                {
                    M(i,i) = Theta[m-1-i];
                    ublas::column(S, m-1-i) = ublas::unit_vector<scalar_type> (S.size1(), i);
                }

                mu = ietl::two_norm(V[0]);
                theta = Theta.back() / mu / mu;
                X_.push_back(V[0]/mu);

                r = W[0] / mu - theta * X_.back();
                norm_r = ietl::two_norm(r);
                
            }//accept
            // restart
            if(m >= iter.m_max()) {
                if(verbose_ > 1) 
                    std::cout<<"restart...\n";
                M = ublas::zero_matrix<scalar_type>(iter.m_min(),iter.m_min());

                vector_set_type tempvs (iter.m_min());
                for(size_t i = 1; i < iter.m_min(); ++i)
                {
                    tempvs[i] = V[0] * S(0,m-1-i);
                    for(size_t j = 1; j < m; ++j)
                        tempvs[i] += V[j] * S(j,m-1-i);
                }
                V.resize(iter.m_min());

                for(size_t i = 1; i < iter.m_min(); ++i)
                    std::swap(V[i],tempvs[i]);

                for(size_t i = 1; i < iter.m_min(); ++i)
                {
                    tempvs[i] = W[0] * S(0,m-1-i);
                    for(size_t j = 1; j < m; ++j)
                        tempvs[i] += W[j] * S(j,m-1-i);
                }
                W.resize(iter.m_min());
                for(size_t i = 1; i < iter.m_min(); ++i)
                    std::swap(W[i],tempvs[i]);

                Theta.resize(iter.m_min());
                for(size_t i = 0; i < iter.m_min(); ++i)
                    M(i,i) = Theta[m-1-i];

                W[0] = w_tilde;
                V[0] = u_tilde;
                m = iter.m_min();

            }// restart

            theta += tau;

            // correction equation
            solver(theta, r, t, iter); //solver is allowed to change r

            //assure t is orthogonal to Q
            norm_w = ietl::two_norm(t);

            for(size_t i = 0; i < X_.size();++i)
                t -= ietl::dot(t,X_[i])*X_[i];

            if(ietl::two_norm(t) <= 0.25 * norm_w)
                for(size_t i = 0; i< X_.size();++i)
                    t -= ietl::dot(t,X_[i])*X_[i];

            ++iter;
            if(iter.error_code() == 1)
                throw std::runtime_error(iter.error_msg());
            if(verbose_ > 1)
                std::cout << "JD iteration " << iter.iterations() << "\tresid = " << norm_r << "\n";         
        } //main loop

    }// jdqr interior exterior

}// end namespace ietl
#endif
