/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN

#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include <boost/mpl/list.hpp>
#include "utils/io.hpp" // has to be first include because of impi
#include <iostream>
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "dmrg/sim/matrix_types.h"
#include "test_mps.h"

typedef boost::mpl::list<
#ifdef HAVE_TwoU1PG
TwoU1PG
#endif
> symmetries;

BOOST_AUTO_TEST_CASE_TEMPLATE( Test_MPS_Compress, S, symmetries )
{
    Fixture<S> f;
    // Test MPS compression
    double norm1 = norm(f.mps);
    BOOST_CHECK_CLOSE(norm1, 1.0, 1e-8);
    f.mps.normalize_left();     // Checks that normalization is preserved
    double norm2 = norm(f.mps);
    BOOST_CHECK_CLOSE(norm2, 1.0, 1e-8);
    // Checks that after compression the normalization is preserved
    double compression_trace;
    auto mps2 = compression::l2r_compress(f.mps, 32, 1e-20, compression_trace);
    mps2.normalize_left();
    double norm3 = norm(mps2);
    BOOST_CHECK_CLOSE(norm3, 1.0, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(compression_trace), 1.0, 1e-8);

    // Check some random values
    BOOST_CHECK_EQUAL(mps2[1].data().n_blocks(), 9);
    BOOST_CHECK_CLOSE(std::abs(mps2[1].data()[4](0,0)), 0.61380232931524104, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[1].data()[4](3,0)), 0.78945975231570409, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[1].data()[4](1,1)), 0.70710678118614889, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[1].data()[4](2,1)), 0.70710678118614845, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[1].data()[4](3,2)), 0.61380232931579326, 1e-8);

    BOOST_CHECK_EQUAL(mps2[2].data().n_blocks(), 14);
    BOOST_CHECK_CLOSE(std::abs(mps2[2].data()[3](0,0)), 0.86358071578894613, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[2].data()[4](0,0)), 0.10708900415704489, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[2].data()[4](2,0)), 0.10693866927266986, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[2].data()[9](0,0)), 0.12106610568661753, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[2].data()[9](3,0)), 0.06822157836692079, 1e-8);

    BOOST_CHECK_EQUAL(mps2[3].data().n_blocks(), 9);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[0](0,0)), 0.3494797405967025, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[0](1,0)), 0.0059432339005499706, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[0](3,0)), 0.001800346775965004, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[1](0,0)), 0.21769837259839703, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[1](11,1)), 0.59209488750654526, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[7](9,1)), 0.45264009928067311, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(mps2[3].data()[8](7,0)), 0.83281243970818863, 1e-8);

    BOOST_CHECK_EQUAL(mps2[4].data().n_blocks(), 4);
}

BOOST_AUTO_TEST_CASE_TEMPLATE( Test_MPS_Rotate, S, symmetries )
{
    Fixture<S> f;
    // Test MPS rotation
    f.mps.canonize(0);
    matrix t = matrix(6,6, {1.0000000000000588, -3.0567539950976905e-14, 2.760512847099443e-07, -2.7058888747042958e-14, -2.0386219974666109e-07, -3.5550050708240602e-13, 3.0567531716797842e-14, 1,
    -1.913451608258381e-14, 2.3600317316650584e-21, 1.448071040659122e-14, 1.0559927431889349e-16, -2.760512719947314e-07, 1.9134523617605794e-14, 0.99999999999996381, -8.5930036151948383e-15, -6.237152668522358e-08,
    -1.0907976004561944e-13, 2.7058886172869296e-14, -2.360032375940198e-21, 8.5930102972685836e-15, 1, 1.2627019978663527e-14, 2.1938496646804564e-20, 2.0386219974664904e-07, -1.4480716638155465e-14, 6.2371582961645389e-08,
    -1.2627026030907235e-14, 0.99999999999997713, 8.5518084530118431e-14, 3.5550050708238511e-13, -1.0559927432976027e-16, 1.0907985818199059e-13, -2.1938507203576664e-20, -8.5518163806719196e-14, 1});

    int scale_inactive = 1;

    mps_rotate::rotate_mps(f.mps, t, scale_inactive);

    // check some random values. Reference values are obtained from the aforementioned OpenMOLCAS calculation
    BOOST_CHECK_CLOSE(f.mps[1].data()[4](0,0), -0.61380217712623453, 1e-8);
    BOOST_CHECK_CLOSE(f.mps[1].data()[4](3,0), -0.78945987064212642, 1e-8);
    BOOST_CHECK_EQUAL(f.mps[1].data().n_blocks(), 9);
    BOOST_CHECK_EQUAL(f.mps[1].data()[3].num_rows(), 2);
    BOOST_CHECK_CLOSE(f.mps[3].data()[0](0,0), -0.34947974240673951, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(f.mps[3].data()[0](3,0)), std::abs(-0.001800346670096781), 1e-8);
    BOOST_CHECK_CLOSE(f.mps[3].data()[3](2,1), -1.4831922467019456e-07, 10);
    BOOST_CHECK_CLOSE(f.mps[3].data()[0](11,0), 0.93692334071306171, 1e-8);
    BOOST_CHECK_CLOSE(std::abs(f.mps[4].data()[2](1,0)), 0.083918122130324874, 1e-4); // TODO: check if sign is correct
    BOOST_CHECK_CLOSE(std::abs(f.mps[4].data()[2](3,0)), 1.0322307057313466e-06, 1e-2);
}
