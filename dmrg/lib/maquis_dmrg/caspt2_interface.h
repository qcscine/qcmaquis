/*
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#ifndef CASPT2_INTERFACE_H
#define CASPT2_INTERFACE_H

#include <iostream>
#include <string>
#include <vector>

#include "maquis_dmrg.h"

namespace maquis {

template <class ValueType> // real or complex
class CASPT2_Interface {
public:
  CASPT2_Interface(const std::vector<double>& epsa);

  void compute_4rdm(const std::vector<ValueType> &epsa);
};

} // namespace maquis

#endif // CASPT2_INTERFACE_H
