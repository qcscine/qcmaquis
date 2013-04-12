
#define BOOST_TEST_MODULE ambient_c_kernels
#include <mpi.h>
#include <boost/test/unit_test.hpp>
#include <boost/test/test_case_template.hpp>
#include <boost/mpl/list.hpp>

#include "alps/numeric/matrix.hpp"
#include "alps/numeric/matrix/algorithms.hpp"
#include "alps/numeric/diagonal_matrix.hpp"
#include "ambient/numeric/matrix.hpp"
#include "utilities.h"

#include <boost/numeric/bindings/lapack/driver/geev.hpp>

template<class T>
struct test_complex_only{
    static void test(int size){}
};

template<class T>
struct test_complex_only<std::complex <T> >{
    static void test(int size){
        typedef alps::numeric::matrix<std::complex<T> > sMatrix;
        typedef ambient::numeric::tiles<ambient::numeric::matrix<std::complex<T> > > pMatrix;
        
        typedef ambient::numeric::tiles<ambient::numeric::diagonal_matrix<std::complex<T> > > pDiagMatrix;
        typedef alps::numeric::diagonal_matrix<std::complex<T> > sDiagMatrix;
        
        pMatrix  pA(size, size);
        pMatrix pLV(size, size);
        pMatrix pRV(size, size);
        
        sMatrix  sA(size, size);
        sMatrix sLV(size, size);
        sMatrix sRV(size, size);
        
        pDiagMatrix pS(size); 
        sDiagMatrix sS(size, 0);
        
        generate(pA);
        sA = maquis::bindings::matrix_cast<sMatrix>(pA);
        
        geev(pA,pLV,pRV,pS);
        
        typename alps::numeric::associated_vector<alps::numeric::matrix< std::complex<T> > >::type Sv(num_rows(sA));
        
        int info = boost::numeric::bindings::lapack::geev('N', 'V', sA, Sv, sLV, sRV);
        if (info != 0)
          throw std::runtime_error("Error in GEEV !");
        
        typename alps::numeric::associated_diagonal_matrix<alps::numeric::matrix< std::complex<T> > >::type S(Sv);
        
        BOOST_CHECK(pS == S);
        BOOST_CHECK(pRV == sRV);
    }
};



BOOST_AUTO_TEST_CASE_TEMPLATE( HEEV_COMPARISON, T, test_types)
{
    test_complex_only<typename T::value_type>::test(T::valuex);
}
