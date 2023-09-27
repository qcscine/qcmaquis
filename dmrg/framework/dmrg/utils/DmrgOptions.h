/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_DMRG_OPTIONS_H
#define UTILS_DMRG_OPTIONS_H

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/parallel/params.hpp"

#include <boost/filesystem/path.hpp>
#include <string>

//=======================================================================
// Options
//
// a class containing the options set by the user
//-----------------------------------------------------------------------

class DmrgOptions
{
public:
    std::string programname;    // name of the executable
    double time_limit;          // time limit for the simulation
    bool valid;                 // shall we really run?
    DmrgParameters parms;       // parameter object
    
    DmrgOptions(int argc, char** argv);
    DmrgOptions();
};

#endif

