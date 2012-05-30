#define BOOST_TEST_MODULE symmetry

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"

typedef U1 SymmGroup;
typedef maquis::types::dense_matrix<double> Matrix;

BOOST_AUTO_TEST_CASE(block_matrix_constructor_one){
    Index<U1> rows,cols;
    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));
    
    block_matrix<Matrix, U1> ba(rows,cols);

    BOOST_CHECK_EQUAL(ba[0](0,0),0);
    BOOST_CHECK_EQUAL(ba[1](0,0),0);
    BOOST_CHECK_EQUAL(true,ba[1](0,0)!=1);
}

BOOST_AUTO_TEST_CASE(block_matrix_constructor_two){
    block_matrix<Matrix, U1> ba(1,1,*(new Matrix(1,1)));
    BOOST_CHECK_EQUAL(ba[0](0,0),0);
}

BOOST_AUTO_TEST_CASE(block_matrix_left_basis_right_basis){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 2));
    rows.insert(std::make_pair(1, 3));
    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 3));
    
    block_matrix<Matrix, U1> ba(rows,cols);
    
    Index<U1> rows_res(ba.left_basis()),cols_res(ba.right_basis());

    BOOST_CHECK_EQUAL(rows, rows_res);
    BOOST_CHECK_EQUAL(cols, cols_res);
}

BOOST_AUTO_TEST_CASE(block_matrix_n_blocks){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 2));
    rows.insert(std::make_pair(1, 3));
    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 3));
    
    block_matrix<Matrix, U1> ba(rows,cols);
    BOOST_CHECK_EQUAL(ba.n_blocks(),3);
}

BOOST_AUTO_TEST_CASE(block_matrix_operators_operator_bracket){
    Matrix m(1,1,1);
    block_matrix<Matrix, U1> ba(1,1,*(new Matrix(1,1)));
    BOOST_CHECK_EQUAL(ba[0](0,0),0);//const []
    ba[0] = m;// non const []
    BOOST_CHECK_EQUAL(ba[0](0,0),1);
}

BOOST_AUTO_TEST_CASE(block_matrix_has_block){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(3, 1));
    rows.insert(std::make_pair(2, 1));
    rows.insert(std::make_pair(1, 1));

    cols.insert(std::make_pair(3, 1));
    cols.insert(std::make_pair(2, 1));
    cols.insert(std::make_pair(1, 1));

    block_matrix<Matrix, U1> ba(rows,cols);

    BOOST_CHECK_EQUAL(true, ba.has_block(3,3));
    BOOST_CHECK_EQUAL(true, ba.has_block(2,2));
    BOOST_CHECK_EQUAL(true, ba.has_block(1,1));
    BOOST_CHECK_EQUAL(false, ba.has_block(1,2));//check a ghost pair
   
    BOOST_CHECK_EQUAL(true, ba.has_block(std::make_pair(1, 1),std::make_pair(1, 1)));  
    BOOST_CHECK_EQUAL(true, ba.has_block(std::make_pair(2, 1),std::make_pair(2, 1)));  
    BOOST_CHECK_EQUAL(true, ba.has_block(std::make_pair(3, 1),std::make_pair(3, 1)));  
    BOOST_CHECK_EQUAL(false, ba.has_block(std::make_pair(2, 1),std::make_pair(3, 1)));  
}

BOOST_AUTO_TEST_CASE(block_matrix_operator_parenthesis){
  // todo 
}


BOOST_AUTO_TEST_CASE(block_matrix_remove_rows_cols_from_block){
    Index<U1> rows,cols;
    rows.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 2));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,0) = 1;
    ba[0](0,1) = 2;
    ba[0](1,0) = 3;
    ba[0](1,1) = 4;

    ba.remove_rows_from_block(0,0,1);
    BOOST_CHECK_EQUAL(ba[0].num_rows(),1);
    ba.remove_cols_from_block(0,1,1);
    BOOST_CHECK_EQUAL(ba[0](0,0),3);
}

BOOST_AUTO_TEST_CASE(block_matrix_remove_trace){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 2));
    rows.insert(std::make_pair(1, 2));
    rows.insert(std::make_pair(1, 2));

    cols.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 2));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,0) = 1;
    ba[0](1,1) = 2;

    ba[1](0,0) = 3;
    ba[1](1,1) = 4;

    ba[2](0,0) = 5;
    ba[2](1,1) = 6;

    BOOST_CHECK_EQUAL(ba.trace(),21);
}

BOOST_AUTO_TEST_CASE(block_matrix_inplace_transpose){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 2));
    rows.insert(std::make_pair(1, 2));

    cols.insert(std::make_pair(1, 2));
    cols.insert(std::make_pair(1, 2));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,1) = 1;
    ba[0](1,0) = 2;

    ba[1](0,1) = 3;
    ba[1](1,0) = 4;
   
    ba.inplace_transpose();

    BOOST_CHECK_EQUAL(ba[0](0,1),2);
    BOOST_CHECK_EQUAL(ba[0](1,0),1);

    BOOST_CHECK_EQUAL(ba[1](0,1),4);
    BOOST_CHECK_EQUAL(ba[1](1,0),3);
}

BOOST_AUTO_TEST_CASE(block_matrix_clear){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 1));

    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));

    block_matrix<Matrix, U1> ba(rows,cols);
    ba.clear();

    BOOST_CHECK_EQUAL(ba.n_blocks(),0);
   //question to michel why not ras.clear() ??
}


BOOST_AUTO_TEST_CASE(block_matrix_operator_multiply_equal){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 1));

    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,0) = 1 ;
    ba[1](0,0) = 2 ;

    ba *= 2;

    BOOST_CHECK_EQUAL(ba[0](0,0),2);
    BOOST_CHECK_EQUAL(ba[1](0,0),4);
}

BOOST_AUTO_TEST_CASE(block_matrix_operator_divide_equal){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 1));
    rows.insert(std::make_pair(1, 1));

    cols.insert(std::make_pair(1, 1));
    cols.insert(std::make_pair(1, 1));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,0) = 4 ;
    ba[1](0,0) = 8 ;

    ba /= 2;

    BOOST_CHECK_EQUAL(ba[0](0,0),2);
    BOOST_CHECK_EQUAL(ba[1](0,0),4);
}

BOOST_AUTO_TEST_CASE(block_matrix_resize_block){
    Index<U1> rows,cols;

    rows.insert(std::make_pair(1, 3));
    rows.insert(std::make_pair(1, 3));

    block_matrix<Matrix, U1> ba(rows,cols);

    ba[0](0,0) = 1 ;
    ba[0](0,1) = 2 ;
    ba[0](0,2) = 3 ;

    ba[0](1,0) = 4 ;
    ba[0](1,1) = 5 ;
    ba[0](1,2) = 6 ;

    ba[0](2,0) = 7 ;
    ba[0](2,1) = 8 ;
    ba[0](2,2) = 9 ;

    ba.resize_block(3,3,2,2,false);
/*

    BOOST_CHECK_EQUAL(ba[0](0,0),1);
    BOOST_CHECK_EQUAL(ba[0](0,1),2);

    BOOST_CHECK_EQUAL(ba[0](1,0),4);
    BOOST_CHECK_EQUAL(ba[0](1,1),5);
*/
}


