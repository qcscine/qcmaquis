/*
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include <iostream>
#include <string>
#include <vector>

#include "caspt2_interface.h"

namespace maquis {

template <class ValueType>
CASPT2_Interface<ValueType>::CASPT2_Interface(const std::vector<double> &epsa) {
  std::cout << "Printing epsa = \n";
  printf("Printing epsa = \n");
  for (const auto &e : epsa) {
    std::cout << e << " ";
  }
  std::cout << "\n" << std::endl;
  exit(1);
}

template class CASPT2_Interface<double>;
} // namespace maquis
