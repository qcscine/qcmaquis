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
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/any.hpp>
#include <utility>

namespace detail {
    class Logger_impl_base
    {
    public:
        virtual void log(boost::any const &) = 0;
        virtual void save(alps::hdf5::archive & ar) const = 0;
    };

    template<class T>
    class Logger_impl : public Logger_impl_base
    {
    public:
        void log(boost::any const & val)
        {
            vals.push_back(boost::any_cast<T>(val));
        }
        
        void save(alps::hdf5::archive & ar) const
        {
            std::vector<T> allvalues;
            if (ar.is_data("mean/value"))
                ar >> alps::make_pvp("mean/value", allvalues);
            std::copy(vals.begin(), vals.end(), std::back_inserter(allvalues));
            ar << alps::make_pvp("mean/value", allvalues);
        }
        
    private:
        std::vector<T> vals;
    };
}

template<class T>
std::pair<std::string, T> make_log(std::string const & name,
                                   T const & val)
{
    return std::make_pair(name, val);
}

class Logger
{
public:
    template<class T>
    void operator<<(std::pair<std::string, T> const & p)
    {
        if (loggers.count(p.first) > 0)
            loggers[p.first]->log(p.second);
        else {
            loggers[p.first].reset(new detail::Logger_impl<T>());
            loggers[p.first]->log(p.second);
        }
    }
    
    void save(alps::hdf5::archive & ar) const
    {
        for (std::map<std::string, boost::shared_ptr<detail::Logger_impl_base> >::const_iterator
             it = loggers.begin(); it != loggers.end(); ++it)
        {
            ar << alps::make_pvp(it->first, *it->second);
        }
    }
    
private:
    std::map<std::string, boost::shared_ptr<detail::Logger_impl_base> > loggers;
};

#endif
