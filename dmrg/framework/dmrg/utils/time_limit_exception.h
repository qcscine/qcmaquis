/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG__UTILS_TIME_LIMIT_EXCEPTION_H
#define MAQUIS_DMRG__UTILS_TIME_LIMIT_EXCEPTION_H

#include <string>
#include <boost/lexical_cast.hpp>
#include <stdexcept>

namespace dmrg {
    class time_limit : public std::runtime_error {
    public:
        time_limit(int sw, int st)
        : sweep_(sw)
        , site_(st)
        , std::runtime_error(  std::string("time limit reached. current status is [ ")
                             + std::string("sweep=") + boost::lexical_cast<std::string>(sweep_) + std::string(", ")
                             + std::string("site=" ) + boost::lexical_cast<std::string>(site_ ) + std::string(" ]."))
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
