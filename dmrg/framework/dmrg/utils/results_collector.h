/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef UTILS_RESULTS_COLLECTOR_H
#define UTILS_RESULTS_COLLECTOR_H

#include <map>
#include <boost/any.hpp>
#include <boost/shared_ptr.hpp>
#include <utility>

#include "dmrg/utils/storage.h"

namespace detail {
    class collector_impl_base
    {
    public:
        virtual void collect(boost::any const &) = 0;
        virtual void save(alps::hdf5::archive & ar) const = 0;
        // TODO: fixed storage type because templated virtual function are not allowed
    };

    template<class T>
    class collector_impl : public collector_impl_base
    {
    public:
        void collect(boost::any const & val)
        {
            vals.push_back(boost::any_cast<T>(val));
        }
        
        void save(alps::hdf5::archive & ar) const
        {
            std::vector<T> allvalues;
            if (ar.is_data("mean/value"))
                ar["mean/value"] >> allvalues;
            std::copy(vals.begin(), vals.end(), std::back_inserter(allvalues));
            ar["mean/value"] << allvalues;
        }
        
    private:
        std::vector<T> vals;
    };
}

class collector_proxy {
    typedef boost::shared_ptr<detail::collector_impl_base> coll_type;
public:
    collector_proxy(coll_type & collector_)
    : collector(collector_)
    { }
    
    template<class T>
    void operator<<(T const& val)
    {
        if (!collector)
            collector.reset(new detail::collector_impl<T>());
        collector->collect(val);
    }
    
private:
    coll_type & collector;
};

class results_collector
{
public:
    
    void clear()
    {
        collection.clear();
    }
    
    collector_proxy operator[] (std::string name)
    {
        return collector_proxy(collection[name]);
    }
    
    template <class Archive>
    void save(Archive & ar) const
    {
        for (std::map<std::string, boost::shared_ptr<detail::collector_impl_base> >::const_iterator
             it = collection.begin(); it != collection.end(); ++it)
        {
            ar[it->first] << *it->second;
        }
    }
    
private:
    std::map<std::string, boost::shared_ptr<detail::collector_impl_base> > collection;
};

#endif
