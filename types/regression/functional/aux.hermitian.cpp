#include "params.hpp"

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

            static void c(matrix<T>& a){
                assert(a.num_rows() == a.num_cols());
                size_t lda = a.num_rows();
                T* ad = current(a);
                T* ar = updated(a);

                for(size_t i = 0; i < a.num_cols(); ++i) 
                    for(size_t j = i+1; j < a.num_rows(); ++j)
                        randomize(ar[i*lda+j]);

                for(size_t i = 0; i < a.num_cols(); ++i) 
                        randomize_diag(ar[i*lda+i]);
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

