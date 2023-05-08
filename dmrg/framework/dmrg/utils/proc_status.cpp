/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */


#include "dmrg/utils/proc_status.h"

#include <stdio.h>
#include <iostream>
#include <fstream>
#include <sstream>

#include <regex>

#if defined(__APPLE__)
#include <mach/task.h>
#include <mach/mach_init.h>
#endif

std::string proc_status_mem() {
#if defined(__linux__) || defined(__linux) || defined(linux) || defined(__gnu_linux__)
    std::ifstream ifs("/proc/self/status");
    if (ifs) {
        std::string res;
        std::regex peak_expr("^VmPeak:	([ ]*)([0-9]+) kB");
        std::regex size_expr("^VmSize:	([ ]*)([0-9]+) kB");
        std::string line;
        while (!ifs.eof()) {
            getline(ifs, line);
            std::smatch what;
            if      (std::regex_match(line, what, peak_expr))
                res += what.str(2) + " ";
            else if (std::regex_match(line, what, size_expr))
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
