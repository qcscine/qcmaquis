/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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


#include "dmrg/utils/proc_status.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <boost/regex.hpp>

#if defined(__APPLE__)
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

std::string proc_status_mem() {
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    std::ifstream ifs("/proc/self/status");
    if (ifs) {
        std::string res;
        boost::regex peak_expr("^VmPeak:	([ ]*)([0-9]+) kB");
        boost::regex size_expr("^VmSize:	([ ]*)([0-9]+) kB");
        std::string line;
        while (!ifs.eof()) {
            getline(ifs, line);
            boost::smatch what;
            if      (boost::regex_match(line, what, peak_expr))
                res += what.str(2) + " ";
            else if (boost::regex_match(line, what, size_expr))
                res += what.str(2) + " ";
        }
        ifs.close();
        return res;
    } else {
        std::cerr << "Cannot open /proc/self/status." << std::endl;
    }
#elif defined(__APPLE__)
    // Inspired by:
    // http://miknight.blogspot.com/2005/11/resident-set-size-in-mac-os-x.html
    struct task_basic_info t_info;
    mach_msg_type_number_t t_info_count = TASK_BASIC_INFO_COUNT;
    task_info(current_task(), TASK_BASIC_INFO, (task_info_t)&t_info, &t_info_count);
    std::stringstream ss;
    ss << t_info.resident_size / 1024;
    ss << " ";
    ss << t_info.virtual_size / 1024;
    return ss.str();
#endif
    return std::string();
}
