#define BOOST_TEST_MODULE general_matrix
#include <boost/test/unit_test.hpp>
#include <boost/mpl/list.hpp>

#include <boost/lambda/lambda.hpp>
#include <boost/ref.hpp>

#include "../general_matrix.hpp"
#include "../matrix_interface.hpp"
#include "../resizable_matrix_interface.hpp"

using namespace blas;


typedef boost::mpl::list<float, double, int, unsigned int, long unsigned int, char> test_types;
// long long unsigned int causes problems in boost::iterator facade

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
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename general_matrix<T>::row_element_iterator, typename general_matrix<T>::row_element_iterator> range(row(a,i));
        T iota(i);
        std::transform(range.first,range.second,range.first,boost::lambda::_1 = boost::lambda::var(iota)++);
    }
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
    for(unsigned int j=0; j<num_columns(a); ++j)
    {
        std::pair<typename general_matrix<T>::column_element_iterator, typename general_matrix<T>::column_element_iterator> range(column(a,j));
        T iota(j);
        std::transform(range.first,range.second,range.first,boost::lambda::_1 = boost::lambda::var(iota)++);
    }
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
    for(unsigned int i=0; i<num_rows(a); ++i)
    {
        std::pair<typename general_matrix<T>::row_element_iterator, typename general_matrix<T>::row_element_iterator> range(row(a,i));
        T iota(i);
        std::transform(range.first,range.second,range.first,boost::lambda::_1 = boost::lambda::var(iota)++);
    }
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
   
    // What happens if an exception is thrown?
    // Remains the matrix unchanged if an exception is thrown during the resize process?
    general_matrix<T> c(a); 
    vector<T> test;
    std::size_t max_size = test.max_size();
    try
    {
        resize(a,max_size,max_size);
    }
    catch(std::bad_alloc& dummy)
    {
        std::cerr<<"Expception"<<std::endl;
    }
    BOOST_CHECK_EQUAL(a,c);
} 
