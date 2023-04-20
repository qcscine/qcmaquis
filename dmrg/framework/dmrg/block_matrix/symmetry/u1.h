/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SYMMETRY_U1_H
#define SYMMETRY_U1_H

#include <iostream>
#include <vector>

#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>


class U1
{
public:

  typedef int charge;
  typedef int subcharge;

  static const charge IdentityCharge = 0;
  static const bool finite = false;

  static charge fuse(charge a, charge b) { return a + b; }
  static charge particleNumber(charge a) { return a; }

  template<int R> static charge fuse(const boost::array<charge, R> &v)
  {
    charge ret = 0;
    for (int i = 0; i < R; i++)
      ret += v[i];
    return ret;
 }
};

template <class Archive>
inline void serialize(Archive & ar, U1::charge & c, const unsigned int version)
{
    ar & c;
}

#endif
