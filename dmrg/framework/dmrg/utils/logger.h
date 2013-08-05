/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef LOGGER_H
#define LOGGER_H

#include <map>
#include <string>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <utility>

template<class Archive>
class Logger
{
public:
    template<class T>
    void operator<<(const T& p)
    {
    }
    
    template<class OtherArchive>
    void save(OtherArchive & ar) const
    {
    }
};

#endif
