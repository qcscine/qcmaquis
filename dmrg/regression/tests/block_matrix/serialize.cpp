#define BOOST_TEST_MODULE block_matrix_serialize

#include "selector.h"

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include <boost/geometry/geometries/adapted/boost_array.hpp>

#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"

#include <iostream>
#include <fstream>
#include <boost/archive/text_iarchive.hpp>
#include <boost/archive/text_oarchive.hpp>

#include <boost/mpl/list.hpp>

//
// List of types SymmGroup for which the block_matrix<Matrix,SymmGroup> is tested
//
typedef boost::mpl::list<
TrivialGroup
, U1
, TwoU1
> test_types;


template <class SymmGroup>
block_matrix<Matrix, SymmGroup> init();

template <>
block_matrix<Matrix, U1> init<U1>()
{
    block_matrix<Matrix, U1> ba;
    ba.insert_block(Matrix(1, 1, 0), 1, 1);
    ba.insert_block(Matrix(2, 2, 1), 0, 0);
    ba.insert_block(Matrix(1, 4, 2), 2, 2);

    return ba;
}

template <>
block_matrix<Matrix, TrivialGroup> init<TrivialGroup>()
{
    typename TrivialGroup::charge C = TrivialGroup::IdentityCharge;
    block_matrix<Matrix, TrivialGroup> ba;
    ba.insert_block(Matrix(5, 6, 4), C, C);
    
    return ba;
}

template <>
block_matrix<Matrix, TwoU1> init<TwoU1>()
{
    typedef typename TwoU1::charge charge;
    charge A, B, C, D;
    A[0]=1; B[1]=1; C[0]=2; D[1]=1; D[0]=3;

    block_matrix<Matrix, TwoU1> ba;
    ba.insert_block(Matrix(1, 1, 0), A, B);
    ba.insert_block(Matrix(2, 2, 1), D, D);
    ba.insert_block(Matrix(1, 4, 2), C, C);
    
    return ba;
}


BOOST_AUTO_TEST_CASE_TEMPLATE( save_load, SymmGroup, test_types )
{
    block_matrix<Matrix, SymmGroup> ba1 = init<SymmGroup>();
    
    std::ofstream ofs("blockmatrix_serialize.dat");
    boost::archive::text_oarchive oa(ofs);
    oa << ba1;
    ofs.close();
    
    block_matrix<Matrix, SymmGroup> ba2;
    std::ifstream ifs("blockmatrix_serialize.dat");
    boost::archive::text_iarchive ia(ifs);
    ia >> ba2;
    ifs.close();

    maquis::cout << ba2;
    
    BOOST_CHECK_EQUAL( ba1.left_basis(),  ba2.left_basis()  );
    BOOST_CHECK_EQUAL( ba1.right_basis(), ba2.right_basis() );
    for (std::size_t block=0; block<ba1.n_blocks(); ++block)
        BOOST_CHECK_EQUAL( ba1[block], ba2[block] );
    
}
