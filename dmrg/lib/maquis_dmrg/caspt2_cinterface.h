/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */
#ifndef CASPT2_CINTERFACE_H
#define CASPT2_CINTERFACE_H

#include "caspt2_interface.h"

extern "C" {
void qcmaquis_interface_caspt2_init(const double *epsa, int nasht);
}

#endif
