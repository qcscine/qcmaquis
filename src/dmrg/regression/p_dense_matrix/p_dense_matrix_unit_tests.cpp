#include "utils/zout.hpp"
#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/p_dense_matrix_algorithms.h"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"

#define BOOST_TEST_MODULE p_dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>
#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#define M_SIZE 1024
#define M_P_SIZE 256
#define NP 2
using namespace blas;

//
// List of types T for which the p_dense_matrix<T> is tested
// (long long unsigned int causes problems in boost::iterator facade)
typedef boost::mpl::list<double> test_types;
typedef ambient::dim3 dim3;

namespace type_pairs
{

struct IntDouble
{
    typedef int first_type;
    typedef double second_type;
};

struct DoubleInt
{
    typedef double first_type;
    typedef int second_type;
};
};

struct AmbientConfig {
    AmbientConfig()   { ambient::init(); }
    ~AmbientConfig()  { ambient::finalize(); }
};


//
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt> test_type_pairs;

template <typename OutputIterator, typename T>
T fill_range_with_numbers(OutputIterator begin, OutputIterator end, T iota)
{
    // Unfortunately we can't use the postincrement operator, due to std:complex<>
    // -> so we have to emulate it's behaviour...
    std::transform(begin,end,begin,boost::lambda::_1 = (boost::lambda::var(iota)+=T(1))-T(1));
    return iota;
}

template <typename T>
T fill_matrix_with_numbers(p_dense_matrix<T>& a)
{
    T iota(0);
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename p_dense_matrix<T>::row_element_iterator, typename p_dense_matrix<T>::row_element_iterator> range(row(a,i));
        iota += fill_range_with_numbers(range.first,range.second,T(i));
    }
    return iota;
}

BOOST_GLOBAL_FIXTURE( AmbientConfig );

/*
BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, T, test_types )
{
    p_dense_matrix<T> a;
    BOOST_CHECK_EQUAL(num_rows(a), 0 );
    BOOST_CHECK_EQUAL(num_columns(a), 0 );
    a.resize(10, 10);
    BOOST_CHECK_EQUAL(num_rows(a), 10 );
    BOOST_CHECK_EQUAL(num_columns(a), 10 );

    p_dense_matrix<T> b(10,10);
    BOOST_CHECK_EQUAL(num_rows(b), 10 );
    BOOST_CHECK_EQUAL(num_columns(b), 10 );
    for(unsigned int i=0; i<10; ++i)
        for(unsigned int j=0; j<10; ++j)
            BOOST_CHECK_EQUAL(b(i,j), T());

    p_dense_matrix<T> c(10,10,5);
    BOOST_CHECK_EQUAL(num_rows(c), 10 );
    BOOST_CHECK_EQUAL(num_columns(c), 10 );
    for(unsigned int i=0; i<10; ++i)
        for(unsigned int j=0; j<10; ++j)
            BOOST_CHECK_EQUAL(c(i,j), T(5));
}*/

BOOST_AUTO_TEST_CASE_TEMPLATE( summ_operation_test, T, test_types )
{
    Timer time1("ambient");
/*    Timer time2("pure");
    time2.begin();

    double** ad = (double**)malloc(sizeof(double*)*M_SIZE/M_P_SIZE/NP);
    for(int i=0; i < M_SIZE/M_P_SIZE/NP; i++){
        ad[i] = (double*)malloc(sizeof(double)*M_P_SIZE*M_P_SIZE);
        memset(ad[i], 0, sizeof(double)*M_P_SIZE*M_P_SIZE);
    }
    double** bd = (double**)malloc(sizeof(double*)*M_SIZE/M_P_SIZE/NP);
    for(int i=0; i < M_SIZE/M_P_SIZE/NP; i++){
        bd[i] = (double*)malloc(sizeof(double)*M_P_SIZE*M_P_SIZE);
        memset(bd[i], 0, sizeof(double)*M_P_SIZE*M_P_SIZE);
    }
    double** cd = (double**)malloc(sizeof(double*)*M_SIZE/M_P_SIZE/NP);
    for(int i=0; i < M_SIZE/M_P_SIZE/NP; i++){
        cd[i] = (double*)malloc(sizeof(double)*M_P_SIZE*M_P_SIZE);
        memset(cd[i], 0, sizeof(double)*M_P_SIZE*M_P_SIZE);
    }

    for(int k=0; k < M_SIZE/M_P_SIZE/NP; k++)
    for(int i=0; i < M_P_SIZE*M_P_SIZE; i++)
    ad[k][i] = bd[k][i] + cd[k][i];
    MPI_Barrier(MPI_COMM_WORLD);
    time2.end();
*/
    ambient::layout >> dim3(10,5), dim3(2,2), dim3(10,1);

    p_dense_matrix<T> a(M_SIZE,M_SIZE);
    p_dense_matrix<T> b(M_SIZE,M_SIZE);
    p_dense_matrix<T> c(M_SIZE,M_SIZE);

    time1.begin();
    a = b * c;
    ambient::push(ambient::block_2d_cyclic_l_kernel, ambient::null_c_kernel, c);
    ambient::playout();

    MPI_Barrier(MPI_COMM_WORLD);
    time1.end();
}

/*BOOST_AUTO_TEST_CASE_TEMPLATE( single_integer_test, T, test_types )
{
    int *input = (int*)malloc(sizeof(int));
    *input = -1;
    ambient::push(ambient::single_integer_l_kernel, ambient::single_integer_c_kernel, *input);
    ambient::playout();
    BOOST_CHECK_EQUAL(*input, 13);
}*/
