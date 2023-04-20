/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
// Substitute header for Boost::UnitTest floating point comparison to avoid warnings

#ifndef FPCOMPARISON_H
#define FPCOMPARISON_H

#ifdef LEGACY_BOOST
#include <boost/test/floating_point_comparison.hpp>
#else
#include <boost/test/tools/floating_point_comparison.hpp>
#endif

#endif