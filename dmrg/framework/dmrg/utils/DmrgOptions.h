/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef UTILS_DMRG_OPTIONS_H
#define UTILS_DMRG_OPTIONS_H

#include "dmrg/utils/DmrgParameters.h"

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

