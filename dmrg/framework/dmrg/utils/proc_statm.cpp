/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include "dmrg/utils/proc_statm.h"

#include <cstdio>
#include <iostream>
#include <fstream>

std::string proc_statm() {
#if defined(__linux__) || defined(__linux) || defined(linux) || \
    defined(__gnu_linux__)
  std::ifstream ifs("/proc/self/statm");
  if (ifs) {
    std::string statm;
    getline(ifs, statm);
    ifs.close();
    return statm;
  } else {
    std::cerr << "Cannot open /proc/self/statm." << std::endl;
  }
#endif
  return std::string();
}
