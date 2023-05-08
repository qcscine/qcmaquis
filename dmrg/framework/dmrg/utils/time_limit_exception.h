/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG__UTILS_TIME_LIMIT_EXCEPTION_H
#define MAQUIS_DMRG__UTILS_TIME_LIMIT_EXCEPTION_H

#include <string>
#include <boost/lexical_cast.hpp>
#include <stdexcept>

namespace dmrg {
    class time_limit : public std::runtime_error {
    public:
        time_limit() : std::runtime_error("time limit reached") {}
        
        time_limit(int sw, int st)
        : std::runtime_error(  std::string("time limit reached. current status is [ ")
                             + std::string("sweep=") + boost::lexical_cast<std::string>(sw) + std::string(", ")
                             + std::string("site=" ) + boost::lexical_cast<std::string>(st) + std::string(" ]."))
        , sweep_(sw)
        , site_(st)
        { }
        
        int sweep() const throw()
        { return sweep_; }

        int site() const throw()
        { return site_; }
        
    private:
        int sweep_, site_;
    };
}

#endif
