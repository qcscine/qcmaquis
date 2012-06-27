#ifndef MT_MATRIX_H
#define MT_MATRIX_H

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"

#ifdef USE_GPU
#include "types/mt_matrix/matrix_gpu_functions.hpp"
#endif

#include <utils/function_objects.h>
#include <utils/timings.h>
#include <utils/data_collector.hpp>

#include "mtm_scheduler.h"

#ifndef GEMM_MT_THRESHOLD
#define GEMM_MT_THRESHOLD 50
#endif
#ifndef GEMM_GPU_THRESHOLD
#define GEMM_GPU_THRESHOLD 100
#endif

namespace alps { 
    namespace numeric {
        
        template<typename T>
        class mt_matrix;
    }
}

namespace detail {
    
    template <class T>
    void gpu_gemm(alps::numeric::matrix<T> const & A,
                  alps::numeric::matrix<T> const & B,
                  alps::numeric::matrix<T> & C)
    {
        if (A.num_cols() > GEMM_GPU_THRESHOLD)
            gpu::matrix_matrix_multiply(A, B, C, 0.);
        else
            gemm(A, B, C);
    }
    
    template <class T>
    void gpu_gemm(alps::numeric::matrix<std::complex<T> > const & A,
                  alps::numeric::matrix<std::complex<T> > const & B,
                  alps::numeric::matrix<std::complex<T> > & C)
    {
        gemm(A, B, C);
    }
    
}

template<typename T>
class MtmGemmRequest : public MtmRequest
{
public:
    MtmGemmRequest(alps::numeric::mt_matrix<T> const & A,
                   alps::numeric::mt_matrix<T> const & B,
                   alps::numeric::mt_matrix<T> & C,
                   boost::shared_ptr<boost::mutex> mutA,
                   boost::shared_ptr<boost::mutex> mutB,
                   boost::shared_ptr<boost::mutex> mutC)
    : A_(A), B_(B) , C_(C)
    , mutA_(mutA), mutB_(mutB), mutC_(mutC)
    , lgA(*mutA), lgB(*mutB), lgC(*mutC)
    { }
    
    void operator()()
    {
        //        A_.wait(); B_.wait(); C_.wait();
        assert( A_.data_.num_cols() == B_.data_.num_rows() );
        static Timer timer("gemm requested");
        timer.begin();
#ifdef USE_GPU
        detail::gpu_gemm(A_.data_, B_.data_, C_.data_);
#else
        gemm(A_.data_, B_.data_, C_.data_);
#endif
        timer.end();
    }
    
private:
    alps::numeric::mt_matrix<T> const & A_;
    alps::numeric::mt_matrix<T> const & B_;
    alps::numeric::mt_matrix<T> & C_;
    boost::shared_ptr<boost::mutex> mutA_, mutB_, mutC_;
    boost::lock_guard<boost::mutex> lgA, lgB, lgC;
};

namespace alps
{
    namespace numeric
    {
        
        template <typename T, class DiagMatrix>
        void svd(mt_matrix<T> const & M, mt_matrix<T> & U, mt_matrix<T> & V,
                             DiagMatrix & S);
        
        template <typename T>
        void qr(mt_matrix<T> const & M, mt_matrix<T> & Q,  mt_matrix<T> & R);
        
        template <typename T>
        void syev(mt_matrix<T> const & M, mt_matrix<T> & vecs,
                              typename alps::numeric::diagonal_matrix<double> & S);

        template <typename T>
        void heev(mt_matrix<T> M, mt_matrix<T> & evecs,
                              typename alps::numeric::associated_real_vector<matrix<T> >::type & evals);
        
        template <typename T>
        void heev(mt_matrix<T> M, mt_matrix<T> & evecs,
                              typename alps::numeric::associated_diagonal_matrix<matrix<T> >::type & evals);
        
        template <typename T>
        mt_matrix<T> exp (mt_matrix<T> M, T const & alpha=1);
        
        
        template<typename T>
        class mt_matrix
        {   
            typedef alps::numeric::matrix<T> slave_t;
            
        public:
            typedef typename slave_t::value_type value_type;
            typedef typename slave_t::reference reference;
            typedef typename slave_t::const_reference const_reference;
            typedef typename slave_t::size_type size_type;
            typedef typename slave_t::difference_type difference_type;
            
            typedef typename slave_t::element_iterator element_iterator;
            
            static mt_matrix<T> identity_matrix (size_type size)
            {
                mt_matrix<T> ret;
                ret.data_ = slave_t::identity_matrix(size);
                return ret;
            }
            
            explicit mt_matrix(size_type rows = 0, size_type cols = 0,
                               value_type val = value_type())
            : data_(rows, cols, val)
            , consistent(true)
            { }
            
            ~mt_matrix()
            {
                wait();
            }
            
            mt_matrix(mt_matrix const & rhs)
            {
                rhs.wait();
                data_ = rhs.data_;
                consistent = true;
            }
            
            void swap_with(mt_matrix & r)
            {
                wait(); r.wait();
                swap(this->data_, r.data_);
            }
            
            friend void swap(mt_matrix & x,
                             mt_matrix & y)
            {
                x.swap_with(y);
            }
            
            mt_matrix& operator=(mt_matrix r)
            {
                this->swap_with(r);
                return *this;
            }
            
            value_type & operator()(size_type i, size_type j)
            {
                wait();
                return data_(i, j);
            }
            
            value_type const & operator()(size_type i, size_type j) const
            {
                wait();
                return data_(i, j);
            }
            
            mt_matrix& operator+=(mt_matrix const & rhs)
            {
                wait(); rhs.wait();
                data_ += rhs.data_;
                return *this;
            }
            
            mt_matrix& operator-=(mt_matrix const & rhs)
            {
                wait(); rhs.wait();
                data_ -= rhs.data_;
                return *this;
            }
            
            template<typename T2>
            mt_matrix& operator*=(T2 t)
            {
                wait();
                data_ *= t;
                return *this;
            }
            
            template<typename T2>
            mt_matrix& operator/=(T2 t)
            {
                wait();
                data_ /= t;
                return *this;
            }
            
            void conjugate_inplace()
            {
                wait();
                conj_inplace(data_);
            }
            
            std::pair<element_iterator,element_iterator> elements()
            {
                wait();
                return data_.elements();
            }
            
            void remove_rows(size_type i, difference_type k)
            {
                wait();
                data_.remove_rows(i, k);
            }
            
            void remove_cols(size_type i, difference_type k)
            {
                wait();
                data_.remove_cols(i, k);
            }
            
            void resize(size_type r, size_type c)
            {
                wait();
                data_.resize(r, c);
            }
            
#ifdef HAVE_ALPS_HDF5
            void load(alps::hdf5::archive & ar)
            {
                wait();
                data_.load(ar);
            }
            void save(alps::hdf5::archive & ar) const
            {
                wait();
                data_.save(ar);
            }
#endif
            
            size_type num_rows() const { wait(); return data_.num_rows(); }
            size_type num_cols() const { wait(); return data_.num_cols(); }
            
            
            friend
            void gemm(mt_matrix const & A,
                      mt_matrix const & B,
                      mt_matrix & C)
            {
                //        gemm(A.data_, B.data_, C.data_);
                
                A.wait(); B.wait(); C.wait();
                if (A.data_.num_cols() < GEMM_MT_THRESHOLD) {
                    static Timer timer("single thread gemm");
                    timer.begin();
                    gemm(A.data_, B.data_, C.data_);
                    timer.end();
                } else {
                    boost::shared_ptr<boost::mutex> mta(new boost::mutex()), mtb(new boost::mutex()), mtc(new boost::mutex());
                    {
                        A.consistent = false;
                        B.consistent = false;
                        C.consistent = false;
                        
                        boost::lock_guard<boost::mutex> lock_listA(A.list_mutex), lock_listB(B.list_mutex), lock_listC(C.list_mutex);
                        A.wait_mutices.push_back(mta);
                        B.wait_mutices.push_back(mtb);
                        C.wait_mutices.push_back(mtc);
                    }
                    
                    DCOLLECTOR_ADD(gemm_collector, A.data_.num_cols())
                    
                    mtm_global_master.push(boost::shared_ptr<MtmRequest>(new MtmGemmRequest<T>(A, B, C,
                                                                                               mta, mtb, mtc)));
                }
                
                //        maquis::cerr << "Waiting..." << std::endl;
                //        A.wait(); B.wait(); C.wait();
                //        maquis::cerr << "Done waiting." << std::endl;
            }
            
            template <typename TT, class DiagMatrix>
            friend
            void algorithms::svd(mt_matrix<TT> const & M, mt_matrix<TT> & U, mt_matrix<TT> & V,
                                 DiagMatrix & S);
            
            template <typename U>
            friend
            void algorithms::qr(mt_matrix<U> const & M, mt_matrix<U> & Q,  mt_matrix<U> & R);
            
            template <typename U>
            friend
            void algorithms::syev(mt_matrix<U> const & M, mt_matrix<U> & vecs,
                                  typename alps::numeric::diagonal_matrix<double> & S);
            template <typename U>
            friend
            void algorithms::heev(mt_matrix<U> M, mt_matrix<U> & evecs,
                      typename alps::numeric::associated_real_vector<matrix<U> >::type & evals);
            
            template <typename U>
            friend
            void algorithms::heev(mt_matrix<U> M, mt_matrix<U> & evecs,
                      typename alps::numeric::associated_diagonal_matrix<matrix<U> >::type & evals);

            template <typename U>
            friend
            mt_matrix<U> algorithms::exp (mt_matrix<U>, U const &);

            
            
            template<class Generator>
            void generate(Generator g)
            {
                wait();
                algorithms::generate(data_, g);
            }
            
            template<typename T2>
            friend mt_matrix operator*(mt_matrix m, T2 t)
            {
                m.wait();
                m *= t;
                return m;
            }
            
            template<typename T2>
            friend mt_matrix operator*(T2 t, mt_matrix m)
            {
                m.wait();
                m *= t;
                return m;
            }
            
        private:
            slave_t data_;
            
            mutable boost::mutex list_mutex;
            mutable std::list<boost::shared_ptr<boost::mutex> > wait_mutices;
            mutable bool consistent;
            
            void wait() const
            {
                if (consistent)
                    return;
                
                boost::lock_guard<boost::mutex> lock_list(list_mutex);
                for (std::list<boost::shared_ptr<boost::mutex> >::iterator it = wait_mutices.begin();
                     it != wait_mutices.end(); ++it)
                    boost::lock_guard<boost::mutex> lock(**it);
                wait_mutices.clear();
                
                consistent = true;
            }
            
            friend class MtmGemmRequest<T>;
        };
        
        template<typename T>
        std::size_t num_rows(mt_matrix<T> const & m)
        {
            return m.num_rows();
        }
        
        template<typename T>
        std::size_t num_cols(mt_matrix<T> const & m)
        {
            return m.num_cols();
        }
        
        template<typename T>
        std::pair<typename mt_matrix<T>::element_iterator, typename mt_matrix<T>::element_iterator>
        elements(mt_matrix<T> & m)
        {
            return m.elements();
        }
        
        template<typename T>
        void resize(mt_matrix<T> & m,
                    typename mt_matrix<T>::size_type r,
                    typename mt_matrix<T>::size_type c)
        {
            m.resize(r, c);
        }
        
        template<typename T>
        mt_matrix<T> transpose(mt_matrix<T> const & m)
        {
            mt_matrix<T> ret(num_cols(m), num_rows(m));
            
            for (std::size_t r = 0; r < num_rows(m); ++r)
                for (std::size_t c = 0; c < num_cols(m); ++c)
                    ret(c, r) = m(r, c);
            return ret;
        }
        
        template<typename T>
        mt_matrix<T> conjugate(mt_matrix<T> m)
        {
            m.conjugate_inplace();
            return m;
        }
        
        template<typename T>
        T trace(mt_matrix<T> const & m)
        {
            T ret = T();
            for (std::size_t r = 0; r < std::min(num_rows(m), num_cols(m)); ++r)
                ret += m(r, r);
            return ret;
        }
        
        template<class T>
        std::size_t size_of(mt_matrix<T> const & m)
        {
            return num_rows(m)*num_cols(m)*sizeof(T);
        }
        
        template<typename T>
        struct alps::numeric::associated_diagonal_matrix<mt_matrix<T> >
        {
            typedef typename alps::numeric::associated_diagonal_matrix<matrix<T> >::type type;
        };
        
        template<typename T>
        struct alps::numeric::associated_vector<mt_matrix<T> >
        {
            typedef typename alps::numeric::associated_vector<matrix<T> >::type type;
        };
        
        template<typename T>
        struct alps::numeric::associated_real_diagonal_matrix<mt_matrix<T> >
        {
            typedef typename alps::numeric::associated_real_diagonal_matrix<matrix<T> >::type type;
        };
        
        template<typename T>
        struct alps::numeric::associated_real_vector<mt_matrix<T> >
        {
            typedef typename alps::numeric::associated_real_vector<matrix<T> >::type type;
        };

    }
}

template<typename T>
std::ostream& operator<<(std::ostream & os, alps::numeric::mt_matrix<T> const & m)
{
    for (std::size_t r = 0; r < num_rows(m); ++r) {
        for (std::size_t c = 0; c < num_cols(m); ++c)
            os << m(r, c);
        os << std::endl;
    }
    return os;
}


#endif
