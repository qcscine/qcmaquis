/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#define BOOST_TEST_MAIN
// Unit test for integral map
#include <boost/test/included/unit_test.hpp>
#ifdef LEGACY_BOOST
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif
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