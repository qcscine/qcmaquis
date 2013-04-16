#include "params.hpp"

template<typename T>
struct random_helper{
    static void randomize(T& a){
        a = drand48();
    }

    static void randomize(std::complex<T>& a){
        a.real(drand48());
        a.imag(drand48());
    }
};

template<class T, std::size_t Size>
struct helper_test{
/*nothing this function is only for TE (complex) */
     static void test(){ }
};

template<class T, std::size_t Size>
struct helper_test<std::complex<T>, Size >{
     static void test(){
         typedef alps::numeric::matrix<std::complex<T>> sMatrix;
         typedef ambient::numeric::tiles<ambient::numeric::matrix<std::complex<T> > > pMatrix;
     
         std::complex<T> a; 
         random_helper<std::complex<T> >::randomize(a);
     
         pMatrix pA(Size, Size), pB(Size, Size);
         sMatrix sA(Size, Size), sB(Size, Size);
        
         generate(pA);
         sA = maquis::bindings::matrix_cast<sMatrix>(pA);
     
         sB = exp(sA,a);
         pB = exp(pA,a);
     
         BOOST_CHECK(pB==sB);
     }
};

BOOST_AUTO_TEST_CASE_TEMPLATE( EXP_ALGO, T, test_types) {
    helper_test<typename T::value_type, T::valuex>::test(); 
}

BOOST_AUTO_TEST_CASE_TEMPLATE( EXP_ALGO_HERMITIAN, T, test_types) {
    typedef alps::numeric::matrix<typename T::value_type> sMatrix;
    typedef ambient::numeric::tiles<ambient::numeric::matrix<typename T::value_type> > pMatrix;

    typename T::value_type a; 
    random_helper<typename T::value_type>::randomize(a);

    int Size = T::valuex; 

    pMatrix pA(Size, Size), pB(Size, Size);
    sMatrix sA(Size, Size), sB(Size, Size);
   
    generate(pA);
    sA = maquis::bindings::matrix_cast<sMatrix>(pA);

    sB = exp_hermitian(sA,a);
    pB = exp_hermitian(pA,a);

    BOOST_CHECK(pB==sB);
}
