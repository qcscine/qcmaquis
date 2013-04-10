
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>
#include <boost/random.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

// C - I am created a hermitian generator maybe to integrate in ambient later
namespace ambient { namespace numeric {

    namespace kernels {

        template<typename T>
        struct init_random_hermitian : public kernel< init_random_hermitian<T> > {

            static void randomize(T& a){
                a = drand48();
            }

            static void randomize(std::complex<T>& a){
                a.real(drand48());
                a.imag(drand48());
            }

            static void randomize_diag(T& a){
                a = drand48();
            }

            static void randomize_diag(std::complex<T>& a){
                a.real(drand48());
                a.imag(0.0);
            }
            static void c(unbound< matrix<T> >& a){
                size_t size = a.num_rows();
                T* ad = updated(a);
                for(size_t i = 0; i < size; ++i) 
                    for(size_t j = i; j < size; ++j){
                        randomize(ad[i*size+j]);
                        ad[j*size+i] = helper_complex<T>::conj(ad[i*size+j]);
                    }

                for(size_t i = 0; i < size; ++i) 
                    randomize_diag(ad[i*size+i]);
            }
        };
    } //end namespace

    template<class Matrix>
    inline void generate_hermitian(tiles<Matrix>& a){
        assert(a.mt == a.nt);
        for(size_t i = 0; i < a.mt; ++i) 
            for(size_t j = i; j < a.nt; ++j){
                fill_random_hermitian(a.tile(i,j));
                a.tile(j,i) = conj(a.tile(i,j));
            }
    }   

    template<class Matrix>
    inline void fill_random_hermitian(Matrix& a){
        kernels::init_random_hermitian<typename Matrix::value_type>::spawn<complexity::N2>(a);
    }   
        
}} //end namespace



BOOST_AUTO_TEST_CASE_TEMPLATE( CONJ_INPLACE, T, test_types)
{
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    pMatrix pA(T::valuex, T::valuey);
    sMatrix sA(T::valuex, T::valuey);

    if(T::valuex == T::valuey)
        generate_hermitian(pA);
    else
        generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    bool bsA = is_hermitian(sA);
    bool bpA = is_hermitian(pA);

    BOOST_CHECK(bpA==bsA);
}

