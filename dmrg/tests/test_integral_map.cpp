/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#define BOOST_TEST_MAIN
// Unit test for integral map
#include <boost/test/included/unit_test.hpp>
#include "utils/fpcomparison.h"
#include <integral_interface.h>

BOOST_AUTO_TEST_CASE( Test_Integral_Map )
{
    // Real integrals with nonrelativistic permutation
    chem::integral_map<double> ints;
    ints[{1,1,2,1}] = 1.0;
    ints[{1,1,1,2}] = 2.0;

    double val = ints[{1,1,2,1}];
    // Two insertions of the symmetry-permuted element must yield only one element in the map
    BOOST_CHECK_EQUAL(ints.size(), 1);

    // check if the permutation works
    BOOST_CHECK_CLOSE(val, 2.0, 1.0e-15);

    // Now try the relativistic permutation
    chem::integral_map<std::complex<double> > ints_complex;
    ints_complex[{1,1,2,1}] = {1.0, 0.0};
    ints_complex[{1,1,1,2}] = {2.0, -1.0};

    // Here the integrals should not permute
    BOOST_CHECK_EQUAL(ints_complex.size(), 2);

    // Check if the insertion actually worked
    BOOST_CHECK_CLOSE(std::real(ints_complex[{1,1,1,2}]), 2.0, 1.0e-15);
    BOOST_CHECK_CLOSE(std::imag(ints_complex[{1,1,1,2}]), -1.0, 1.0e-15);

}