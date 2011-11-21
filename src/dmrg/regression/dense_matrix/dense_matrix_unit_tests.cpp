#define BOOST_TEST_MODULE maquis::types::dense_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <boost/lambda/lambda.hpp>
#include <complex>
#include <numeric>

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"




//
// List of types T for which the maquis::types::dense_matrix<T> is tested
//
typedef boost::mpl::list<float, double, int, unsigned int, long unsigned int,std::complex<float>, std::complex<double> > test_types;
// long long unsigned int causes problems in boost::iterator facade


namespace type_pairs
{
struct DComplexDouble
{
    typedef std::complex<double> first_type;
    typedef double second_type;
};

struct DoubleDComplex
{
    typedef double first_type;
    typedef std::complex<double> second_type;
};

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

//
// List of type pairs <T,U> for which the mixed type matrix vector multiplication is tested.
//
typedef boost::mpl::list<type_pairs::IntDouble, type_pairs::DoubleInt, type_pairs::DoubleDComplex, type_pairs::DComplexDouble> test_type_pairs;

template <typename OutputIterator, typename T>
T fill_range_with_numbers(OutputIterator begin, OutputIterator end, T iota)
{
    // Unfortunately we can't use the postincrement operator, due to std:complex<>
    // -> so we have to emulate it's behaviour...
    std::transform(begin,end,begin,boost::lambda::_1 = (boost::lambda::var(iota)+=T(1))-T(1));
    return iota;
}

template <typename T>
T fill_matrix_with_numbers(maquis::types::dense_matrix<T>& a)
{
    T iota(0);
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename maquis::types::dense_matrix<T>::row_element_iterator, typename maquis::types::dense_matrix<T>::row_element_iterator> range(row(a,i));
        iota += fill_range_with_numbers(range.first,range.second,T(i));
    }
    return iota;
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, T, test_types )
{
    maquis::types::dense_matrix<T> a;
    BOOST_CHECK_EQUAL(num_rows(a), 0 );
    BOOST_CHECK_EQUAL(num_cols(a), 0 );

    maquis::types::dense_matrix<T> b(10,10);
    BOOST_CHECK_EQUAL(num_rows(b), 10 );
    BOOST_CHECK_EQUAL(num_cols(b), 10 );
    for(unsigned int i=0; i<10; ++i)
        for(unsigned int j=0; j<10; ++j)
            BOOST_CHECK_EQUAL(b(i,j), T());

    maquis::types::dense_matrix<T> c(15,5,5);
    BOOST_CHECK_EQUAL(num_rows(c), 15 );
    BOOST_CHECK_EQUAL(num_cols(c), 5 );
    for(unsigned int i=0; i<15; ++i)
        for(unsigned int j=0; j<5; ++j)
            BOOST_CHECK_EQUAL(c(i,j), T(5));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_swap_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(10,10,1);
    maquis::types::dense_matrix<T> b(1,1,0);
    maquis::types::dense_matrix<T> c(a);
    maquis::types::dense_matrix<T> d(b);
    std::swap(a,b);
    BOOST_CHECK_EQUAL(a,d);
    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignement_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(10,10,1);
    maquis::types::dense_matrix<T> b(1,1,0);
    b = a;
    BOOST_CHECK_EQUAL(a,b);
    b(0,0) = 5;
    BOOST_CHECK_EQUAL(a(0,0) != b(0,0), true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( row_iterator_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(10,20);
    fill_matrix_with_numbers(a);

    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename maquis::types::dense_matrix<T>::row_element_iterator, typename maquis::types::dense_matrix<T>::row_element_iterator> range(row(a,i));
        unsigned int j=0;
        for(typename maquis::types::dense_matrix<T>::const_row_element_iterator it(range.first); it != range.second; ++it)
        {
            BOOST_CHECK_EQUAL(a(i,j), *it);
            ++j;
        }
    }
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
            BOOST_CHECK_EQUAL(a(i,j),T(i+j));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( column_iterator_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(10,20);
    fill_matrix_with_numbers(a);
    for(unsigned int j=0; j<num_cols(a); ++j)
    {
        std::pair<typename maquis::types::dense_matrix<T>::column_element_iterator, typename maquis::types::dense_matrix<T>::column_element_iterator> range(column(a,j));
        unsigned int i=0;
        for(typename maquis::types::dense_matrix<T>::const_column_element_iterator it(range.first); it != range.second; ++it)
        {
            BOOST_CHECK_EQUAL(a(i,j), *it);
            ++i;
        }
    }
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
            BOOST_CHECK_EQUAL(a(i,j),T(i+j));
}

BOOST_AUTO_TEST_CASE_TEMPLATE( element_iterator_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(10,20);
    maquis::types::dense_matrix<T> b(10,20);
    std::pair<typename maquis::types::dense_matrix<T>::element_iterator,typename maquis::types::dense_matrix<T>::element_iterator> range(elements(a));
    fill_range_with_numbers(range.first,range.second,0);

    T k = T(0);
    T sum = T(0);
    for(unsigned int j=0; j<num_cols(a); ++j)
        for(unsigned int i=0; i<num_rows(a); ++i)
        {
            b(i,j) = k;
            sum += k;
            k += T(1);
        }

    T acc = std::accumulate(range.first, range.second,T(0));
    BOOST_CHECK_EQUAL(acc,sum);
    BOOST_CHECK_EQUAL(a,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( resize_test, T, test_types )
{
    maquis::types::dense_matrix<T> a;

    // Check primitive enlargement
    resize(a,10,5);
    BOOST_CHECK_EQUAL(num_rows(a),10);
    BOOST_CHECK_EQUAL(num_cols(a),5);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    // Resize case 1:
    // Enlargement out of the reserved range
    // size1 > reserved_size1_
    unsigned int size1 = a.capacity().first + 10;
    // Check whether enlargement keeps the values of the original matrix
    resize(a,size1,15,1);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( i >=10 || j >= 5)
                BOOST_CHECK_EQUAL(a(i,j),T(1));
            else
                BOOST_CHECK_EQUAL(a(i,j),T(i+j));
        }

    // Resize case 2:
    // Shrinking
    // size1 < reserved_size1
    // size1 < size1_ (-> shrinking)
    resize(a,10,5);
    BOOST_CHECK_EQUAL(a,b);

    // Resize case 3:
    // Enlargement within the already reserved range
    // size1 < reserved_size1
    // size1 > size1_
    resize(a,15,10);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( i >= 10 || j >= 5) BOOST_CHECK_EQUAL(a(i,j),T(0));
            else BOOST_CHECK_EQUAL(a(i,j), T(i+j));
        }
   

}

BOOST_AUTO_TEST_CASE_TEMPLATE( resize_exception_test, T, test_types )
{
    maquis::types::dense_matrix<T> a(22,18);
    fill_matrix_with_numbers(a);
    
    // What happens if an exception is thrown?
    // Remains the matrix unchanged if an exception is thrown during the resize process?
    // Case 1: size1 > reserved_size1_
    maquis::types::dense_matrix<T> ref(a); 
    maquis::types::dense_matrix<T> c(a); 
    maquis::types::dense_matrix<T> d(a); 
    vector<T> test;
    std::size_t max_size = test.max_size();
    try
    {
        resize(a,max_size+10,1);
    }
    catch(...)
    {
        BOOST_CHECK_EQUAL(a,ref);
    }

    // Resize case 2:
    // Shrinking in one dimension
    // size1 < reserved_size1
    // size1 < size1_ (-> shrinking)
    try
    {
        resize(c,1,max_size+10);
    }
    catch(...)
    {
        BOOST_CHECK_EQUAL(c,ref);
    }

    // Resize case 3:
    // Enlargement within the already reserved range
    // size1 < reserved_size1
    // size1 > size1_
    resize(d,2,5);
    maquis::types::dense_matrix<T> ref_d(d);
    try
    {
        resize(d,4,max_size/2+5);
    }
    catch(...)
    {
        BOOST_CHECK_EQUAL(d,ref_d);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( reserve_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(22,18);
    fill_matrix_with_numbers(a);
    
    maquis::types::dense_matrix<T> ref(a);

    // Case 1:
    // size1 > reserved_size1_
    a.reserve(30,30);
    BOOST_CHECK_EQUAL(a,ref);
    BOOST_CHECK_EQUAL(a.capacity().first >= 30 && a.capacity().second >= 30, true);
    
    // Case 2:
    // size1 < reserved_size1_
    // reserved_size1_*size2 > values_.capacity

    a.reserve(20,40);
    BOOST_CHECK_EQUAL(a,ref);
    BOOST_CHECK_EQUAL(a.capacity().first >= 30 && a.capacity().second >= 40, true);


    // Case 3:
    // size1 < reserved_size1_
    // size2 < size2_
    a.reserve(10,10);
    BOOST_CHECK_EQUAL(a,ref);
    BOOST_CHECK_EQUAL(a.capacity().first >= 30 && a.capacity().second >= 40, true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( append_rows_test, T, test_types)
{
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    vector<T> data_single(initsize,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    iota = fill_range_with_numbers(data_single.begin(),data_single.end(),iota);
    iota = fill_range_with_numbers(data_multiple.begin(),data_multiple.end(),iota);

    // Append a single row
    append_rows(a, std::make_pair(data_single.begin(), data_single.end()) );
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( i != initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),T(j));
        }
    // Append multiple rows
    append_rows(a, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( i < initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
            {
                switch (i)
                {
                    case initsize:
                        BOOST_CHECK_EQUAL(a(i,j),T(j));
                        break;
                    case initsize+1:
                        BOOST_CHECK_EQUAL(a(i,j),T(j+initsize));
                        break;
                    case initsize+2:
                        BOOST_CHECK_EQUAL(a(i,j),T(j+2*initsize));
                        break;
                    case initsize+3:
                        BOOST_CHECK_EQUAL(a(i,j),T(j+3*initsize));
                        break;
                    default:
                        // There should not be any other row
                        // Report an error
                        BOOST_CHECK( true == false);
                }
            }

        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( append_columns_test, T, test_types)
{
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    vector<T> data_single(initsize,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    iota = fill_range_with_numbers(data_single.begin(),data_single.end(),iota);
    iota = fill_range_with_numbers(data_multiple.begin(),data_multiple.end(),iota);

    // Append a single column
    append_columns(a, std::make_pair(data_single.begin(), data_single.end()) );
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( j != initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),T(i));
        }
    // Append multiple rows
    append_columns(a, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if( j < initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
            {
                switch (j)
                {
                    case initsize:
                        BOOST_CHECK_EQUAL(a(i,j),T(i));
                        break;
                    case initsize+1:
                        BOOST_CHECK_EQUAL(a(i,j),T(i+initsize));
                        break;
                    case initsize+2:
                        BOOST_CHECK_EQUAL(a(i,j),T(i+2*initsize));
                        break;
                    case initsize+3:
                        BOOST_CHECK_EQUAL(a(i,j),T(i+3*initsize));
                        break;
                    default:
                        // There should not be any other column
                        // Report an error
                        BOOST_CHECK( true == false);
                }
            }
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_rows_test, T, test_types)
{
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    // remove the last row
    remove_rows(a,initsize-1);
    // remove the first row
    remove_rows(a,0);
    //remove some rows in the middle
    remove_rows(a,5);
    remove_rows(a,11,4);

    BOOST_CHECK_EQUAL(num_rows(a),initsize-7);

    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if(i<5)
                BOOST_CHECK_EQUAL(a(i,j),b(i+1,j));
            else if (i < 11)
                BOOST_CHECK_EQUAL(a(i,j),b(i+2,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),b(i+6,j));
        }
    
    maquis::types::dense_matrix<T> c(b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_columns_test, T, test_types)
{
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    // remove the last row
    remove_columns(a,initsize-1);
    // remove the first row
    remove_columns(a,0);
    //remove some columns in the middle
    remove_columns(a,5);
    remove_columns(a,11,4);

    BOOST_CHECK_EQUAL(num_cols(a),initsize-7);

    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            if(j<5)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+1));
            else if (j < 11)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+2));
            else
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+6));
        }
    
    maquis::types::dense_matrix<T> c(b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( insert_rows_test, T, test_types)
{
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    vector<T> data_single(20,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    iota = fill_range_with_numbers(data_single.begin(),data_single.end(),iota);
    iota = fill_range_with_numbers(data_multiple.begin(),data_multiple.end(),iota);

    // Insert a row in for the 0th line, the last line and in the middle
    insert_rows(a, initsize, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 0, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 5, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 8, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            switch(i)
            {
                case 0:
                case 5:
                case 25:
                    BOOST_CHECK_EQUAL(a(i,j),T(j));
                    break;
                case 8:
                case 9:
                case 10:
                    BOOST_CHECK_EQUAL(a(i,j),T(j+(i-7)*initsize));
                    break;
                default:
                    if( i>10 )
                        BOOST_CHECK_EQUAL(a(i,j),b(i-5,j));
                    else if( i>5 )
                        BOOST_CHECK_EQUAL(a(i,j),b(i-2,j));
                    else
                        BOOST_CHECK_EQUAL(a(i,j),b(i-1,j));
            }
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( insert_columns_test, T, test_types)
{ 
    const unsigned int initsize = 20;
    maquis::types::dense_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    vector<T> data_single(20,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    iota = fill_range_with_numbers(data_single.begin(),data_single.end(),iota);
    iota = fill_range_with_numbers(data_multiple.begin(),data_multiple.end(),iota);
    
    // Insert a column in for the 0th line, the last line and in the middle
    insert_columns(a, initsize, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 0, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 5, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 8, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
        {
            switch(j)
            {
                case 0:
                case 5:
                case 25:
                    BOOST_CHECK_EQUAL(a(i,j),T(i));
                    break;
                case 8:
                case 9:
                case 10:
                    BOOST_CHECK_EQUAL(a(i,j),T(i+(j-7)*initsize));
                    break;
                default:
                    if( j>10 )
                        BOOST_CHECK_EQUAL(a(i,j),b(i,j-5));
                    else if( j>5 )
                        BOOST_CHECK_EQUAL(a(i,j),b(i,j-2));
                    else
                        BOOST_CHECK_EQUAL(a(i,j),b(i,j-1));
            }
        }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_assign_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);
    
    a += b;
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), T((i+j)*2) );

    a += a;
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), T((i+j)*4) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    maquis::types::dense_matrix<T> zero(20,30,T(0));
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);
    a += b;
    a -= b;
    BOOST_CHECK_EQUAL(a,b);
    
    a -= a;
    BOOST_CHECK_EQUAL(a,zero);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    maquis::types::dense_matrix<T> zero(20,30,T(0));
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);
    a *= T(1);
    BOOST_CHECK_EQUAL(a,b);
    a *= T(0);
    BOOST_CHECK_EQUAL(a,zero);

    fill_matrix_with_numbers(a);
    a *= T(2);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_cols(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), T(i+j)*T(2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);

    maquis::types::dense_matrix<T> c = a + b;
    a +=b;
    BOOST_CHECK_EQUAL(c,a);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);
    a += b;
    maquis::types::dense_matrix<T> c = a - b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    maquis::types::dense_matrix<T> b(a);
    maquis::types::dense_matrix<T> ref_b(b);
    a*= T(2);
    maquis::types::dense_matrix<T> c = T(2) * b;
    //TODO Do we really want to assume commutative types?
    maquis::types::dense_matrix<T> d = b * T(2);
    BOOST_CHECK_EQUAL(c,a);
    BOOST_CHECK_EQUAL(d,a);
    BOOST_CHECK_EQUAL(b,ref_b);

    // Check whether or not it works with mixed types.
    // (value_type != T2 ) - at least for non integer types...
    maquis::types::dense_matrix<T> e(b);
    b*= 5;
    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            typename maquis::types::dense_matrix<T>::value_type tmp (e(i,j));
            tmp *= 5;
            BOOST_CHECK_EQUAL(b(i,j),tmp);
        }
    maquis::types::dense_matrix<T> ref_e(e);
    maquis::types::dense_matrix<T> f ( e * 5 );
    maquis::types::dense_matrix<T> g ( 5 * e );
    BOOST_CHECK_EQUAL(b,f);
    BOOST_CHECK_EQUAL(b,g);
    BOOST_CHECK_EQUAL(ref_e,e);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_vector_multiply_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    vector<T> v(30);
    fill_matrix_with_numbers(a);
    fill_range_with_numbers(v.begin(),v.end(),T(0));
    maquis::types::dense_matrix<T> a_(a);
    vector<T> v_(v);
    
    vector<T> result(a*v);
    BOOST_CHECK_EQUAL(result.size(),num_rows(a));
    BOOST_CHECK_EQUAL(a,a_);
    BOOST_CHECK_EQUAL(v,v_);
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        T row_result(0);
        for(unsigned int j=0; j<num_cols(a); ++j)
            row_result += a(i,j)*v(j);
        BOOST_CHECK_EQUAL(result(i),row_result);
    }

}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_vector_multiply_mixed_types_test, TPair, test_type_pairs)
{
    // -maquis::types::dense_matrix<T> * vector<int>
    
    maquis::types::dense_matrix<typename TPair::first_type> a(20,30);
    vector<typename TPair::second_type> v(30);
    fill_matrix_with_numbers(a);
    fill_range_with_numbers(v.begin(),v.end(),0);
    maquis::types::dense_matrix<typename TPair::first_type> a_(a);
    vector<typename TPair::second_type> v_(v);
    
    vector<typename maquis::types::MultiplyReturnType<typename TPair::first_type,std::vector<typename TPair::first_type>,typename TPair::second_type, std::vector<typename TPair::second_type> >::value_type> result(a*v);
    BOOST_CHECK_EQUAL(result.size(),num_rows(a));
    BOOST_CHECK_EQUAL(a,a_);
    BOOST_CHECK_EQUAL(v,v_);
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        typename maquis::types::MultiplyReturnType<typename TPair::first_type, std::vector<typename TPair::first_type>,typename TPair::second_type, std::vector<typename TPair::second_type> >::value_type row_result(0);
        for(unsigned int j=0; j<num_cols(a); ++j)
            row_result += a(i,j)*v(j);
        BOOST_CHECK_EQUAL(result(i),row_result);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_matrix_multiply_test, T, test_types)
{
    maquis::types::dense_matrix<T> a(20,30);
    maquis::types::dense_matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    maquis::types::dense_matrix<T> c = a * b;

    BOOST_CHECK_EQUAL(num_rows(c), num_rows(a));
    BOOST_CHECK_EQUAL(num_cols(c), num_cols(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_cols(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_cols(a); ++k)
                result += a(i,k) * b(k,j);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}


