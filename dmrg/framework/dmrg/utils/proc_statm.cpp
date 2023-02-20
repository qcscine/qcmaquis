/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */


#include "dmrg/utils/proc_statm.h"

#include <stdio.h>
#include <iostream>
#include <fstream>

std::string proc_statm() {
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
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
