/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2020 by Alberto Baiardi <abaiardi@ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
 *
 *****************************************************************************/

#define BOOST_TEST_MODULE block_matrix

#include <boost/test/included/unit_test.hpp>
#include <boost/mpl/assert.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/block_matrix/block_matrix.h"

#include "dmrg/sim/matrix_types.h"

//typedef alps::numeric::associated_real_diagonal_matrix<matrix>::type DiagMatrix;

#if defined(HAVE_TwoU1) || defined(HAVE_TwoU1PG)

#ifdef MAQUIS_OPENMP
/* add test for omp */
BOOST_AUTO_TEST_CASE(BlockMatrixCheckOMP) {
    block_matrix<matrix, TwoU1> ba;

    #pragma omp parallel for schedule(dynamic)
    for (int i = 0; i < 100; ++i) {
        matrix m(1, 1, i*1.);
        typename TwoU1::charge chargeL(i);
        typename TwoU1::charge chargeR(i+1);
        #pragma omp critical
        {
            ba.insert_block(m, chargeL, chargeR);
        }
    }
    //#pragma omp critical

    BOOST_CHECK_EQUAL(ba.n_blocks(), 100);
}
#endif

/* Checks that the number of overall blocks of a block_matrix is correct */
BOOST_AUTO_TEST_CASE(BlockMatrixCheckConstructorFromIndexesViaBlocks) {
    Index<TwoU1> rows, cols;
    typename TwoU1::charge charge1(0);
    typename TwoU1::charge charge2(1);
    typename TwoU1::charge charge3(2);
    rows.insert(std::make_pair(charge1, 1));
    rows.insert(std::make_pair(charge2, 1));
    rows.insert(std::make_pair(charge3, 1));
    cols.insert(std::make_pair(charge1, 1));
    cols.insert(std::make_pair(charge2, 1));
    cols.insert(std::make_pair(charge3, 1));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    BOOST_CHECK_EQUAL(ba.n_blocks(), 3);
}

/* Checks that with the block_matrix initialization works correctly */
BOOST_AUTO_TEST_CASE(BlockMatrixConstructorFromIndexesViaValues) {
    Index<TwoU1> rows, cols;
    typename TwoU1::charge leftCharge(0);
    typename TwoU1::charge rightCharge(1);
    rows.insert(std::make_pair(leftCharge, 1));
    rows.insert(std::make_pair(rightCharge, 1));
    cols.insert(std::make_pair(leftCharge, 1));
    cols.insert(std::make_pair(rightCharge, 1));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    BOOST_CHECK_EQUAL(ba[0](0,0), 0);
    BOOST_CHECK_EQUAL(ba[0](0,0), 0);
    BOOST_CHECK_EQUAL(true, ba[1](0,0) != 1.);
}

/* Check that the initializer from the DualIndex object works correctly */
BOOST_AUTO_TEST_CASE(BlockMatrixConstructorFromDualIndex) {
    Index<TwoU1> rows, cols;
    typename TwoU1::charge firstCharge(0);
    typename TwoU1::charge secondCharge(1);
    auto block1 = typename dual_index_detail::QnBlock<TwoU1>(firstCharge, firstCharge, 1, 1);
    auto block2 = typename dual_index_detail::QnBlock<TwoU1>(secondCharge, secondCharge, 1, 1);
    DualIndex<TwoU1> dualIndex;
    dualIndex.insert(block1);
    block_matrix<matrix, TwoU1> ba(dualIndex);
    BOOST_CHECK_EQUAL(ba.n_blocks(), 1);
    dualIndex.insert(block2);
    ba = block_matrix<matrix, TwoU1>(dualIndex);
    BOOST_CHECK_EQUAL(ba.n_blocks(), 2);
}

/* Checks that the getter for the index is working */
BOOST_AUTO_TEST_CASE(BlockMatrixCheckIndexGetter) {
    Index<TwoU1> rows, cols;
    typename TwoU1::charge firstCharge(0);
    typename TwoU1::charge secondCharge(1);
    typename TwoU1::charge thirdCharge(3);
    rows.insert(std::make_pair(firstCharge, 1));
    rows.insert(std::make_pair(secondCharge, 5));
    rows.insert(std::make_pair(thirdCharge, 4));
    cols.insert(std::make_pair(firstCharge, 1));
    cols.insert(std::make_pair(secondCharge, 4));
    cols.insert(std::make_pair(thirdCharge, 5));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    auto rows_res = ba.left_basis();
    auto cols_res = ba.right_basis();
    BOOST_CHECK_EQUAL(rows, rows_res);
    BOOST_CHECK_EQUAL(cols, cols_res);
}

/* Checks that the number of elements is calculated correctly */
BOOST_AUTO_TEST_CASE(BlockMatrixCheckNumerOfElements) {
    Index<TwoU1> rows, cols;
    typename TwoU1::charge firstCharge(0);
    typename TwoU1::charge secondCharge(1);
    rows.insert(std::make_pair(firstCharge, 10));
    rows.insert(std::make_pair(secondCharge, 5));
    cols.insert(std::make_pair(firstCharge, 4));
    cols.insert(std::make_pair(secondCharge, 6));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    auto numberOfElements = ba.num_elements();
    BOOST_CHECK_EQUAL(numberOfElements, 70);
}

/* Check on the operator bracket */
BOOST_AUTO_TEST_CASE(BlockMatrixOperatorBracket) {
    Index<TwoU1> rows, cols;
    matrix m(1, 1, 31.);
    rows.insert(std::make_pair(typename TwoU1::charge(0), 1));
    cols.insert(std::make_pair(typename TwoU1::charge(0), 1));
    block_matrix<matrix, TwoU1> ba(rows,cols);
    BOOST_CHECK_EQUAL(ba[0](0,0), 0.);
    ba[0] = m;
    BOOST_CHECK_EQUAL(ba[0](0,0), 31);
}

/* Check on the scaling of a block_matrix times a scalar */
BOOST_AUTO_TEST_CASE(BlockMatrixScaleByScalar) {
    Index<TwoU1> rows,cols;
    rows.insert(std::make_pair(typename TwoU1::charge(0), 1));
    rows.insert(std::make_pair(typename TwoU1::charge(1), 3));
    cols.insert(std::make_pair(typename TwoU1::charge(0), 1));
    cols.insert(std::make_pair(typename TwoU1::charge(1), 3));
    block_matrix<matrix, TwoU1> ba(rows,cols);
    ba[0](2,2) = 1;
    ba[1](0,0) = 2;
    ba *= 2;
    BOOST_CHECK_EQUAL(ba[0](2, 2), 2);
    BOOST_CHECK_EQUAL(ba[1](0, 0), 4);
}

/* Checks division by a scalar */
BOOST_AUTO_TEST_CASE(BlockMatrixDivideByScalar) {
    Index<TwoU1> rows,cols;
    rows.insert(std::make_pair(typename TwoU1::charge(1), 1));
    rows.insert(std::make_pair(typename TwoU1::charge(0), 3));
    cols.insert(std::make_pair(typename TwoU1::charge(1), 1));
    cols.insert(std::make_pair(typename TwoU1::charge(0), 3));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0, 0) = 2.;
    ba[1](0, 0) = 2.;
    ba[1](1, 1) = 4.;
    ba[1](2, 2) = 8.;
    ba /= 2.;
    BOOST_CHECK_EQUAL(ba[0](0,0), 1.);
    BOOST_CHECK_EQUAL(ba[1](1,1), 2.);
    BOOST_CHECK_EQUAL(ba[1](2,2), 4.);
}

/* Checks the has_block function */
BOOST_AUTO_TEST_CASE(BlockMatrixHasBlock){
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(0);
    auto charge1 = typename TwoU1::charge(1);
    auto charge2 = typename TwoU1::charge(2);
    // Note that here we check the insert method with pos (to force un-sorting)
    rows.insert(std::make_pair(charge0, 3));
    rows.insert(std::make_pair(charge1, 2));
    rows.insert(std::make_pair(charge2, 1));
    cols.insert(std::make_pair(charge2, 3));
    cols.insert(std::make_pair(charge1, 2));
    cols.insert(std::make_pair(charge0, 1));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    BOOST_CHECK_EQUAL(true, ba.has_block(charge0, charge0));
    BOOST_CHECK_EQUAL(true, ba.has_block(charge1, charge1));
    BOOST_CHECK_EQUAL(true, ba.has_block(std::make_pair(charge2, 1), std::make_pair(charge2, 3)));
    BOOST_CHECK_EQUAL(true, ba.has_block(std::make_pair(charge2, 3), std::make_pair(charge2, 1)));  
    BOOST_CHECK_EQUAL(false, ba.has_block(charge1, charge2));
}

/* Addition of two block_matrix objects */
BOOST_AUTO_TEST_CASE(BlockMatrixAddition) {
    Index<TwoU1> rows0, rows1, cols0, cols1;
    auto charge0 = typename TwoU1::charge(0);
    auto charge1 = typename TwoU1::charge(1);
    rows0.insert(std::make_pair(charge0, 5));
    rows1.insert(std::make_pair(charge1, 1));
    cols0.insert(std::make_pair(charge0, 1));
    cols1.insert(std::make_pair(charge1, 5));
    block_matrix<matrix, TwoU1> ba0(rows0, cols0), ba1(rows1,cols1);
    //std::cout << "before add: ba1 is " << ba1 << std::endl;
    ba1 += ba0;
    //std::cout << "after  add: ba1 is " << ba1 << std::endl;
    BOOST_CHECK_EQUAL(true, ba1.has_block(charge0, charge0));
}

/* Subtraction of two block_matrix objects */
BOOST_AUTO_TEST_CASE(BlockMatrixSubtraction) {
    Index<TwoU1> rows0, rows1, cols0, cols1;
    auto charge0 = typename TwoU1::charge(0);
    auto charge1 = typename TwoU1::charge(1);
    rows0.insert(std::make_pair(charge0, 5));
    rows0.insert(std::make_pair(charge1, 1));
    rows1.insert(std::make_pair(charge1, 1));
    cols0.insert(std::make_pair(charge0, 1));
    cols0.insert(std::make_pair(charge1, 5));
    cols1.insert(std::make_pair(charge1, 5));
    block_matrix<matrix, TwoU1> ba0(rows0, cols0), ba1(rows1,cols1);
    ba0[0](0, 0) = 10.;
    ba0[0](0, 1) = 20.;
    ba0[0](0, 4) = 30.;
    ba0[1](0, 0) = 100.;

    ba1[0](0, 0) = 30.;
    ba1[0](0, 1) = 20.;
    ba1[0](0, 4) = 10.;
    //std::cout << "before sub: ba0 is " << ba0 << std::endl;
    ba0 -= ba1;
    //std::cout << "after  sub: ba0 is " << ba0 << std::endl;
    //std::cout << "after  sub: ba1 is " << ba1 << std::endl;
    BOOST_CHECK_EQUAL(ba0[0](0,0), -20.);
    BOOST_CHECK_EQUAL(ba0[0](0,1), 0.);
    BOOST_CHECK_EQUAL(ba0[0](0,4), 20.);
    BOOST_CHECK_EQUAL(ba0[1](0,0), 100.);
}

/* Shift the basis of the block_matrix */
BOOST_AUTO_TEST_CASE(BlockMatrixChargeShift) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(4);
    auto charge1 = typename TwoU1::charge(3);
    auto charge2 = typename TwoU1::charge(2);
    auto charge3 = typename TwoU1::charge(1);
    rows.insert(std::make_pair(charge0, 5));
    rows.insert(std::make_pair(charge1, 4));
    rows.insert(std::make_pair(charge2, 3));
    rows.insert(std::make_pair(charge3, 2));
    cols.insert(std::make_pair(charge0, 1));
    cols.insert(std::make_pair(charge1, 2));
    cols.insert(std::make_pair(charge2, 3));
    cols.insert(std::make_pair(charge3, 4));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 1;
    //std::cout << "BEFORE SHIFT: ba is " << ba << std::endl;
    ba.shift_basis(charge3);
    ba.shift_basis(-charge2);
    //std::cout << "AFTER  SHIFT: ba is " << ba << std::endl;
    BOOST_CHECK_EQUAL(false, ba.has_block(charge1, charge2));
    BOOST_CHECK_EQUAL(false, ba.has_block(charge0, charge0));
    BOOST_CHECK_EQUAL(true, ba.has_block(charge1, charge1));
}


/* Test on the find_block method */
BOOST_AUTO_TEST_CASE(BlockMatrixFindBlock) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(4);
    auto charge1 = typename TwoU1::charge(3);
    auto charge2 = typename TwoU1::charge(2);
    auto charge3 = typename TwoU1::charge(1);
    auto charge4 = typename TwoU1::charge(0);
    rows.insert(std::make_pair(charge0, 5));
    rows.insert(std::make_pair(charge1, 4));
    rows.insert(std::make_pair(charge2, 3));
    rows.insert(std::make_pair(charge3, 2));
    rows.insert(std::make_pair(charge4, 1));
    cols.insert(std::make_pair(charge0, 1));
    cols.insert(std::make_pair(charge1, 2));
    cols.insert(std::make_pair(charge2, 3));
    cols.insert(std::make_pair(charge3, 4));
    cols.insert(std::make_pair(charge4, 5));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    //std::cout << "ba is " << ba << std::endl;
    //std::cout << "find c3 c3 is " << charge3 << " <-> " << charge3 << " == " << ba.find_block(charge3,charge3) << std::endl;
    //std::cout << "find c3 c4 is " << charge0 << " <-> " << charge0 << " == " << ba.find_block(charge0,charge0) << std::endl;
    //std::cout << "find c3 c4 is " << charge3 << " <-> " << charge4 << " == " << ba.find_block(charge3,charge4) << std::endl;
    //std::cout << "find c3 c4 is " << charge0 << " <-> " << charge4 << " == " << ba.find_block(charge0,charge4) << std::endl;
    /* FIND_BLOCK returns the data.size() AKA a "pointer" to the last elemenet for pairs of charges (charge3,charge4) that are 
    not in the block matrix...*/
    BOOST_CHECK_EQUAL(ba.find_block(charge3,charge3), 3);
    BOOST_CHECK_EQUAL(ba.find_block(charge0,charge0), 0);
    BOOST_CHECK_EQUAL(ba.find_block(charge3,charge4), 5);
    BOOST_CHECK_EQUAL(ba.find_block(charge0,charge4), 5);
}

/* Test on the insert_block method */
BOOST_AUTO_TEST_CASE(BlockMatrixInsertBlock) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(4);
    auto charge1 = typename TwoU1::charge(3);
    auto charge2 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 5));
    rows.insert(std::make_pair(charge1, 1));
    cols.insert(std::make_pair(charge1, 5));
    cols.insert(std::make_pair(charge0, 1));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    //std::cout << "ba is " << ba << std::endl;
    matrix m(1, 1, 10.);
    ba.insert_block(m, charge0, charge2);
    //std::cout << "ba is " << ba << std::endl;
    BOOST_CHECK_EQUAL(ba[1](0,0), 10.);
}

/* Test method to remove a block from the block_matrix object */
BOOST_AUTO_TEST_CASE(BlockMatrixRemoveBlock){
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(3);
    auto charge1 = typename TwoU1::charge(2);
    auto charge2 = typename TwoU1::charge(1);
    auto charge3 = typename TwoU1::charge(0);
    rows.insert(std::make_pair(charge0, 30));
    rows.insert(std::make_pair(charge1, 20));
    rows.insert(std::make_pair(charge2, 10));
    rows.insert(std::make_pair(charge3, 1));
    cols.insert(std::make_pair(charge0, 30));
    cols.insert(std::make_pair(charge1, 20));
    cols.insert(std::make_pair(charge2, 10));
    cols.insert(std::make_pair(charge3, 1));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 1.;
    ba[1](0,0) = 2.;
    ba[2](0,0) = 3.;
    ba[3](0,0) = 4.;
    ba.remove_block(1);
    BOOST_CHECK_EQUAL(ba[0](0,0), 1.);
    BOOST_CHECK_EQUAL(ba[1](0,0), 3.);
    ba.remove_block(charge3, charge3);
    BOOST_CHECK_EQUAL(ba[0](0,0), 1.);
    BOOST_CHECK_EQUAL(ba[1](0,0), 3.);
    ba.remove_block(charge0, charge0);
    BOOST_CHECK_EQUAL(ba[0](0,0), 3.);
}

/* Checks that the calculation of the trace is correct */
BOOST_AUTO_TEST_CASE(BlockMatrixCalculateTrace){
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(0);
    auto charge1 = typename TwoU1::charge(1);
    auto charge2 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 6));
    rows.insert(std::make_pair(charge1, 5));
    rows.insert(std::make_pair(charge2, 7));
    cols.insert(std::make_pair(charge0, 6));
    cols.insert(std::make_pair(charge1, 5));
    cols.insert(std::make_pair(charge2, 7));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 1.;
    ba[0](1,1) = 2.;
    ba[0](1,2) = -4.;
    ba[1](0,0) = 3.;
    ba[1](1,1) = 6.;
    ba[1](1,4) = 2.;
    ba[2](0,0) = 5.;
    ba[2](1,1) = 6.;
    BOOST_CHECK_EQUAL(ba.trace(), 23);
}

/* Check on the norm method, which interprets the block_matrix as a 1D vector */
BOOST_AUTO_TEST_CASE(BlockMatrixCalculateNorm) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(0);
    auto charge2 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 2));
    rows.insert(std::make_pair(charge2, 3));
    cols.insert(std::make_pair(charge0, 3));
    cols.insert(std::make_pair(charge2, 4));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 3.;
    ba[0](1,1) = 4.;
    ba[0](2,1) = 10.;
    ba[1](0,1) = 10.;
    //std::cout << "ba is " << ba << std::endl;
    //std::cout << "norm is " << ba.norm() << std::endl;
    BOOST_CHECK_EQUAL(ba.norm(), 15);
}

/* Checks block_matrix transposition */
BOOST_AUTO_TEST_CASE(BlockMatrixTransposeInPlace) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(0);
    auto charge2 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 1));
    rows.insert(std::make_pair(charge2, 1));
    cols.insert(std::make_pair(charge0, 3));
    cols.insert(std::make_pair(charge2, 3));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 1.;
    ba[0](0,1) = 2.;
    ba[0](0,2) = 3.;
    ba[1](0,0) = 4.;
    ba[1](0,1) = 5.;
    ba[1](0,2) = 6.;
    ba.transpose_inplace();
    BOOST_CHECK_EQUAL(ba.left_basis()[0].second, 3);
    BOOST_CHECK_EQUAL(ba.right_basis()[1].second, 1);
    BOOST_CHECK_EQUAL(ba[0](1,0), 2.);
    BOOST_CHECK_EQUAL(ba[1](2,0), 6.);
}

/* Setting to zero of a block_matrix */
BOOST_AUTO_TEST_CASE(BlockMatrixClear) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(0);
    auto charge2 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 1));
    rows.insert(std::make_pair(charge2, 1));
    cols.insert(std::make_pair(charge0, 3));
    cols.insert(std::make_pair(charge2, 3));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba.clear();
    BOOST_CHECK_EQUAL(ba.n_blocks(), 0);
}

/**
 * Test the method to add a block to the block_matrix 
 * not fitting the previous dimensions 
 */
BOOST_AUTO_TEST_CASE(BlockMatrixMatchAndAddBlock) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(3);
    auto charge1 = typename TwoU1::charge(2);
    auto charge2 = typename TwoU1::charge(1);
    auto charge3 = typename TwoU1::charge(0);
    rows.insert(std::make_pair(charge0, 3));
    rows.insert(std::make_pair(charge1, 2));
    rows.insert(std::make_pair(charge2, 1));
    rows.insert(std::make_pair(charge3, 4));
    cols.insert(std::make_pair(charge0, 3));
    cols.insert(std::make_pair(charge1, 2));
    cols.insert(std::make_pair(charge2, 1));
    cols.insert(std::make_pair(charge3, 4));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0, 0) = 1.;
    ba[0](1, 1) = 1.;
    ba[0](2, 2) = 1.;
    ba[0](0, 1) = 2.;
    ba[0](0, 2) = 2.;
    ba[0](1, 2) = 2.;
    ba[0](1, 0) = 2.;
    ba[0](2, 0) = 2.;
    ba[0](2, 1) = 2.;
    ba[1](0, 0) = 3.;
    ba[1](1, 1) = 3.;
    ba[1](0, 1) = 4.;
    ba[1](1, 0) = 4.;
    ba[2](0, 0) = 10.;
    ba[3](0, 0) = 100.;
    ba[3](1, 1) = 100.;
    ba[3](2, 2) = 100.;
    ba[3](3, 3) = 100.;
    // Add matrix to the first block --> smaller matrix
    matrix M1(2, 2, 4.);
    ba.match_and_add_block(M1, charge0, charge0);
    BOOST_CHECK_EQUAL(ba[0](0, 0), 5.);
    BOOST_CHECK_EQUAL(ba[0](0, 1), 6.);
    BOOST_CHECK_EQUAL(ba[0](0, 2), 2.);
    BOOST_CHECK_EQUAL(ba[0](2, 2), 1.);
    // Add matrix to the second block --> same size
    matrix M2(2, 2, -1.);
    ba.match_and_add_block(M2, charge1, charge1);
    BOOST_CHECK_EQUAL(ba[1](0, 0), 2.);
    BOOST_CHECK_EQUAL(ba[1](1, 1), 2.);
    // Add matrix to the third block --> larger
    matrix M3(3, 3, 10.);
    ba.match_and_add_block(M3, charge2, charge2);
    BOOST_CHECK_EQUAL(ba[2](0, 0), 20.);
    BOOST_CHECK_EQUAL(ba[2](0, 1), 10.);
    BOOST_CHECK_EQUAL(ba[2](2, 1), 10.);
    BOOST_CHECK_EQUAL(ba[2](2, 2), 10.);
}

/* Test on method to remove elements that are smaller than a given threshold */
BOOST_AUTO_TEST_CASE(BlockMatrixCleanupZeros) {
    // Remove a zero block in the middle
    auto charge0 = typename TwoU1::charge(3);
    auto charge1 = typename TwoU1::charge(2);
    auto charge2 = typename TwoU1::charge(1);
    block_matrix<matrix, TwoU1> ba;
    matrix m1(4, 6, 1.);
    m1(3,2) = m1(2,2) = m1(3,3) = 1e-15;
    ba.insert_block(m1, charge0, charge0);
    matrix m2(3, 2, 1e-15);
    ba.insert_block(m2, charge1, charge1);
    matrix m3(2, 1, 1.);
    m3(1,0) = 1e-15;
    ba.insert_block(m3, charge2, charge2);
    ba.cleanup_zeros(1e-10);
    BOOST_CHECK_EQUAL(ba.n_blocks(), 2);
    BOOST_CHECK_EQUAL(ba(charge0, charge0)(3,2), 0.);
    BOOST_CHECK_EQUAL(ba(charge0, charge0)(2,2), 0.);
    BOOST_CHECK_EQUAL(ba(charge0, charge0)(3,3), 0.);
    BOOST_CHECK_EQUAL(ba(charge2, charge2)(1,0), 0.);
    ba.cleanup_zeros(10.);
    BOOST_CHECK_EQUAL(ba.n_blocks(),0);
}

/* Reserve blocks in the block_matrix (not used I think) */
BOOST_AUTO_TEST_CASE(BlockMatrixReserve) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(3);
    auto charge1 = typename TwoU1::charge(2);
    auto charge2 = typename TwoU1::charge(1);
    rows.insert(std::make_pair(charge0, 2));
    cols.insert(std::make_pair(charge0, 2));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    BOOST_CHECK_EQUAL(true, ba.has_block(charge0, charge0));
    // Reserve new blocks
    ba.reserve(charge1, charge1, 3, 4);
    BOOST_CHECK_EQUAL(ba.left_basis()[1].second, 3);
    BOOST_CHECK_EQUAL(ba.right_basis()[1].second, 4);
    ba.reserve(charge2, charge2 ,2, 3);
    BOOST_CHECK_EQUAL(true, ba.has_block(charge2, charge2));
    BOOST_CHECK_EQUAL(ba.left_basis()[2].second, 2);
    BOOST_CHECK_EQUAL(ba.left_basis()[2].first, charge2);
    // Reserve old blocks (== resize)
    ba.reserve(charge0, charge0, 3, 1);
    BOOST_CHECK_EQUAL(ba.left_basis()[0].second, 3);
    BOOST_CHECK_EQUAL(ba.right_basis()[0].second, 2);
}

/* Test the allocation routine (to be used together with reserve) */
BOOST_AUTO_TEST_CASE(BlockMatrixAllocate) {
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(3);
    rows.insert(std::make_pair(charge0, 2));
    cols.insert(std::make_pair(charge0, 2));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    ba[0](0,0) = 1.;
    ba[0](0,1) = 1.;
    ba[0](1,0) = 1.;
    ba[0](1,1) = 1.;
    ba.reserve(charge0, charge0, 4, 4);
    ba.allocate_blocks();
    BOOST_CHECK_EQUAL(ba[0](1,1), 1.);
    BOOST_CHECK_EQUAL(ba[0](2,2), 0.);
}

/* Test method to resize a block */
BOOST_AUTO_TEST_CASE(BlockMatrixResizeBlock){
    Index<TwoU1> rows, cols;
    auto charge0 = typename TwoU1::charge(3);
    auto charge1 = typename TwoU1::charge(2);
    rows.insert(std::make_pair(charge0, 10));
    rows.insert(std::make_pair(charge1, 5));
    cols.insert(std::make_pair(charge0, 30));
    cols.insert(std::make_pair(charge1, 3));
    block_matrix<matrix, TwoU1> ba(rows, cols);
    for (int iRow = 0; iRow < 10; iRow++)
        for (int iCol = 0; iCol < 30; iCol++)
            ba[0](iRow, iCol) = 1.;
    ba.resize_block(charge0 , charge0, 20, 30, false);
    BOOST_CHECK_EQUAL(ba[0](0, 0), 1.);
    BOOST_CHECK_EQUAL(ba[0](0, 19), 1.);
    BOOST_CHECK_EQUAL(ba[0](0, 29), 1.);
    BOOST_CHECK_EQUAL(ba[0](9, 19), 1.);
    BOOST_CHECK_EQUAL(ba[0](19, 29), 0.);
}
#endif

