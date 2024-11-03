/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include "caspt2_cinterface.h"

std::unique_ptr<maquis::CASPT2_Interface<double>> caspt2_interface_ptr;

extern "C" {
typedef double V;
void qcmaquis_caspt2_init(const double* epsa, int nasht) {
  std::vector<double> epsa_vec(epsa, epsa + nasht);
  caspt2_interface_ptr.reset(
      new maquis::CASPT2_Interface<double>(epsa_vec));
}
}
