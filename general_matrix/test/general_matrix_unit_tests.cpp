#define BOOST_TEST_MODULE general_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <boost/lambda/lambda.hpp>

#include "../general_matrix.hpp"
#include "../general_matrix_blas.hpp"
#include "../matrix_interface.hpp"
#include "../resizable_matrix_interface.hpp"

using namespace blas;


typedef boost::mpl::list<float, double, int, unsigned int, long unsigned int> test_types;
// long long unsigned int causes problems in boost::iterator facade

template <typename T>
void fill_matrix_with_numbers(general_matrix<T>& a)
{
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename general_matrix<T>::row_element_iterator, typename general_matrix<T>::row_element_iterator> range(row(a,i));
        T iota(i);
        std::transform(range.first,range.second,range.first,boost::lambda::_1 = boost::lambda::var(iota)++);
    }
}

BOOST_AUTO_TEST_CASE_TEMPLATE( constructors_test, T, test_types )
{
    general_matrix<T> a;
    BOOST_CHECK_EQUAL(num_rows(a), 0 );
    BOOST_CHECK_EQUAL(num_columns(a), 0 );

    general_matrix<T> b(10,10);
    BOOST_CHECK_EQUAL(num_rows(b), 10 );
    BOOST_CHECK_EQUAL(num_columns(b), 10 );
    for(unsigned int i=0; i<10; ++i)
        for(unsigned int j=0; j<10; ++j)
            BOOST_CHECK_EQUAL(b(i,j), T());

    general_matrix<T> c(15,5,5);
    BOOST_CHECK_EQUAL(num_rows(c), 15 );
    BOOST_CHECK_EQUAL(num_columns(c), 5 );
    for(unsigned int i=0; i<15; ++i)
        for(unsigned int j=0; j<5; ++j)
            BOOST_CHECK_EQUAL(c(i,j), 5);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( copy_swap_test, T, test_types )
{
    general_matrix<T> a(10,10,1);
    general_matrix<T> b(1,1,0);
    general_matrix<T> c(a);
    general_matrix<T> d(b);
    std::swap(a,b);
    BOOST_CHECK_EQUAL(a,d);
    BOOST_CHECK_EQUAL(b,c);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( assignement_test, T, test_types )
{
    general_matrix<T> a(10,10,1);
    general_matrix<T> b(1,1,0);
    b = a;
    BOOST_CHECK_EQUAL(a,b);
    b(0,0) = 5;
    BOOST_CHECK_EQUAL(a(0,0) != b(0,0), true);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( row_iterator_test, T, test_types )
{
    general_matrix<T> a(10,20);
    fill_matrix_with_numbers(a);

    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename general_matrix<T>::row_element_iterator, typename general_matrix<T>::row_element_iterator> range(row(a,i));
        unsigned int j=0;
        for(typename general_matrix<T>::const_row_element_iterator it(range.first); it != range.second; ++it)
        {
            BOOST_CHECK_EQUAL(a(i,j), *it);
            ++j;
        }
    }
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
            BOOST_CHECK_EQUAL(a(i,j),i+j);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( column_iterator_test, T, test_types )
{
    general_matrix<T> a(10,20);
    fill_matrix_with_numbers(a);
    for(unsigned int j=0; j<num_columns(a); ++j)
    {
        std::pair<typename general_matrix<T>::column_element_iterator, typename general_matrix<T>::column_element_iterator> range(column(a,j));
        unsigned int i=0;
        for(typename general_matrix<T>::const_column_element_iterator it(range.first); it != range.second; ++it)
        {
            BOOST_CHECK_EQUAL(a(i,j), *it);
            ++i;
        }
    }
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
            BOOST_CHECK_EQUAL(a(i,j),i+j);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( resize_test, T, test_types )
{
    general_matrix<T> a;

    // Check primitive enlargement
    resize(a,10,5);
    BOOST_CHECK_EQUAL(num_rows(a),10);
    BOOST_CHECK_EQUAL(num_columns(a),5);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    // Resize case 1:
    // Enlargement out of the reserved range
    // size1 > reserved_size1_
    unsigned int size1 = a.capacity().first + 10;
    // Check whether enlargement keeps the values of the original matrix
    resize(a,size1,15,1);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( i >=10 || j >= 5)
                BOOST_CHECK_EQUAL(a(i,j),1);
            else
                BOOST_CHECK_EQUAL(a(i,j),i+j);
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
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( i >= 10 || j >= 5) BOOST_CHECK_EQUAL(a(i,j),0);
            else BOOST_CHECK_EQUAL(a(i,j), i+j);
        }
   

}

BOOST_AUTO_TEST_CASE_TEMPLATE( resize_exception_test, T, test_types )
{
    general_matrix<T> a(22,18);
    fill_matrix_with_numbers(a);
    
    // What happens if an exception is thrown?
    // Remains the matrix unchanged if an exception is thrown during the resize process?
    // Case 1: size1 > reserved_size1_
    general_matrix<T> ref(a); 
    general_matrix<T> c(a); 
    general_matrix<T> d(a); 
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
    general_matrix<T> ref_d(d);
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
    general_matrix<T> a(22,18);
    fill_matrix_with_numbers(a);
    
    general_matrix<T> ref(a);

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
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    vector<T> data_single(initsize,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    std::transform(data_single.begin(),data_single.end(),data_single.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);
    std::transform(data_multiple.begin(),data_multiple.end(),data_multiple.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);

    // Append a single row
    append_rows(a, std::make_pair(data_single.begin(), data_single.end()) );
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( i != initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),j);
        }
    // Append multiple rows
    append_rows(a, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( i < initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
            {
                switch (i)
                {
                    case initsize:
                        BOOST_CHECK_EQUAL(a(i,j),j);
                        break;
                    case initsize+1:
                        BOOST_CHECK_EQUAL(a(i,j),j+initsize);
                        break;
                    case initsize+2:
                        BOOST_CHECK_EQUAL(a(i,j),j+2*initsize);
                        break;
                    case initsize+3:
                        BOOST_CHECK_EQUAL(a(i,j),j+3*initsize);
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
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    vector<T> data_single(initsize,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    std::transform(data_single.begin(),data_single.end(),data_single.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);
    std::transform(data_multiple.begin(),data_multiple.end(),data_multiple.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);

    // Append a single column
    append_columns(a, std::make_pair(data_single.begin(), data_single.end()) );
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( j != initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),i);
        }
    // Append multiple rows
    append_columns(a, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if( j < initsize)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j));
            else
            {
                switch (j)
                {
                    case initsize:
                        BOOST_CHECK_EQUAL(a(i,j),i);
                        break;
                    case initsize+1:
                        BOOST_CHECK_EQUAL(a(i,j),i+initsize);
                        break;
                    case initsize+2:
                        BOOST_CHECK_EQUAL(a(i,j),i+2*initsize);
                        break;
                    case initsize+3:
                        BOOST_CHECK_EQUAL(a(i,j),i+3*initsize);
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
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    // remove the last row
    remove_rows(a,initsize-1);
    // remove the first row
    remove_rows(a,0);
    //remove some rows in the middle
    remove_rows(a,5);
    remove_rows(a,11,4);

    BOOST_CHECK_EQUAL(num_rows(a),initsize-7);

    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if(i<5)
                BOOST_CHECK_EQUAL(a(i,j),b(i+1,j));
            else if (i < 11)
                BOOST_CHECK_EQUAL(a(i,j),b(i+2,j));
            else
                BOOST_CHECK_EQUAL(a(i,j),b(i+6,j));
        }
    
    general_matrix<T> c(b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( remove_columns_test, T, test_types)
{
    const unsigned int initsize = 20;
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    // remove the last row
    remove_columns(a,initsize-1);
    // remove the first row
    remove_columns(a,0);
    //remove some columns in the middle
    remove_columns(a,5);
    remove_columns(a,11,4);

    BOOST_CHECK_EQUAL(num_columns(a),initsize-7);

    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            if(j<5)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+1));
            else if (j < 11)
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+2));
            else
                BOOST_CHECK_EQUAL(a(i,j),b(i,j+6));
        }
    
    general_matrix<T> c(b);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( insert_rows_test, T, test_types)
{
    const unsigned int initsize = 20;
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    vector<T> data_single(20,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    std::transform(data_single.begin(),data_single.end(),data_single.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);
    std::transform(data_multiple.begin(),data_multiple.end(),data_multiple.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);

    // Insert a row in for the 0th line, the last line and in the middle
    insert_rows(a, initsize, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 0, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 5, std::make_pair(data_single.begin(), data_single.end()) );
    insert_rows(a, 8, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            switch(i)
            {
                case 0:
                case 5:
                case 25:
                    BOOST_CHECK_EQUAL(a(i,j),j);
                    break;
                case 8:
                case 9:
                case 10:
                    BOOST_CHECK_EQUAL(a(i,j),j+(i-7)*initsize);
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
    general_matrix<T> a(initsize,initsize);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    vector<T> data_single(20,1);
    vector<T> data_multiple(3*initsize,2);
    T iota(0);
    std::transform(data_single.begin(),data_single.end(),data_single.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);
    std::transform(data_multiple.begin(),data_multiple.end(),data_multiple.begin(),boost::lambda::_1 = boost::lambda::var(iota)++);
    
    // Insert a column in for the 0th line, the last line and in the middle
    insert_columns(a, initsize, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 0, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 5, std::make_pair(data_single.begin(), data_single.end()) );
    insert_columns(a, 8, std::make_pair(data_multiple.begin(),data_multiple.end()),3);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
        {
            switch(j)
            {
                case 0:
                case 5:
                case 25:
                    BOOST_CHECK_EQUAL(a(i,j),i);
                    break;
                case 8:
                case 9:
                case 10:
                    BOOST_CHECK_EQUAL(a(i,j),i+(j-7)*initsize);
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
    general_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);
    
    a += b;
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), (i+j)*2 );

    a += a;
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), (i+j)*4 );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_assign_test, T, test_types)
{
    general_matrix<T> a(20,30);
    general_matrix<T> zero(20,30,T(0));
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);
    a += b;
    a -= b;
    BOOST_CHECK_EQUAL(a,b);
    
    a -= a;
    BOOST_CHECK_EQUAL(a,zero);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_assign_test, T, test_types)
{
    general_matrix<T> a(20,30);
    general_matrix<T> zero(20,30,T(0));
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);
    a *= T(1);
    BOOST_CHECK_EQUAL(a,b);
    a *= T(0);
    BOOST_CHECK_EQUAL(a,zero);

    fill_matrix_with_numbers(a);
    a *= T(2);
    for(unsigned int i=0; i<num_rows(a); ++i)
        for(unsigned int j=0; j<num_columns(a); ++j)
            BOOST_CHECK_EQUAL( a(i,j), T(i+j)*T(2) );
}

BOOST_AUTO_TEST_CASE_TEMPLATE( plus_test, T, test_types)
{
    general_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);

    general_matrix<T> c = a + b;
    a +=b;
    BOOST_CHECK_EQUAL(c,a);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( minus_test, T, test_types)
{
    general_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);
    a += b;
    general_matrix<T> c = a - b;
    BOOST_CHECK_EQUAL(c,b);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( multiplies_test, T, test_types)
{
    general_matrix<T> a(20,30);
    fill_matrix_with_numbers(a);
    general_matrix<T> b(a);
    general_matrix<T> ref_b(b);
    a*= T(2);
    general_matrix<T> c = T(2) * b;
    //TODO Do we really want to assume commutative types?
    general_matrix<T> d = b * T(2);
    BOOST_CHECK_EQUAL(c,a);
    BOOST_CHECK_EQUAL(d,a);
    BOOST_CHECK_EQUAL(b,ref_b);

    // Check whether or not it works with mixed types.
    // (value_type != T2 ) - at least for non integer types...
    general_matrix<T> e(b);
    b*= 5;
    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_columns(c); ++j)
        {
            typename general_matrix<T>::value_type tmp (e(i,j));
            tmp *= 5;
            BOOST_CHECK_EQUAL(b(i,j),tmp);
        }
    general_matrix<T> ref_e(e);
    general_matrix<T> f ( e * 5 );
    general_matrix<T> g ( 5 * e );
    BOOST_CHECK_EQUAL(b,f);
    BOOST_CHECK_EQUAL(b,g);
    BOOST_CHECK_EQUAL(ref_e,e);

}

BOOST_AUTO_TEST_CASE_TEMPLATE( matrix_matrix_multiply_test, T, test_types)
{
    general_matrix<T> a(20,30);
    general_matrix<T> b(30,50);
    fill_matrix_with_numbers(a);
    fill_matrix_with_numbers(b);

    general_matrix<T> c = a * b;

    BOOST_CHECK_EQUAL(num_rows(c), num_rows(a));
    BOOST_CHECK_EQUAL(num_columns(c), num_columns(b));

    for(unsigned int i=0; i<num_rows(c); ++i)
        for(unsigned int j=0; j<num_columns(c); ++j)
        {
            T result(0);
            for(unsigned int k=0; k< num_columns(a); ++k)
                result += a(i,k) * b(k,j);
            BOOST_CHECK_EQUAL(c(i,j),result);
        }
}


