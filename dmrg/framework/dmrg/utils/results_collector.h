/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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
        virtual ~collector_impl_base() {}
        virtual void collect(boost::any const &) = 0;
        virtual void save(alps::hdf5::archive & ar) const = 0;
        virtual void load(alps::hdf5::archive & ar) = 0;
        virtual const std::vector<boost::any>& get() const = 0;
        // TODO: fixed storage type because templated virtual function are not allowed
    };

    template<class T>
    class collector_impl : public collector_impl_base
    {
    public:
        void collect(boost::any const & val)
        {
            vals.push_back(val);
        }

        void save(alps::hdf5::archive & ar) const
        {
            std::vector<T> allvalues;
            if (ar.is_data("mean/value"))
                ar["mean/value"] >> allvalues;
            allvalues.reserve(allvalues.size()+vals.size());
            for(auto&& val : vals)
                allvalues.push_back(boost::any_cast<T>(val));
            ar["mean/value"] << allvalues;
        }

        void load(alps::hdf5::archive & ar)
        {
            // overwrite the current vector
            vals.clear();
            // read from file
            std::vector<T> allvalues;
            if (ar.is_data("mean/value"))
                ar["mean/value"] >> allvalues;
            vals.reserve(allvalues.size());
            for(auto&& val : allvalues)
                vals.push_back(val);
        }

        // TODO: Copying is inefficient!
        const std::vector<boost::any>& get() const { return vals; };

    private:
        std::vector<boost::any> vals;
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

    // Needed for a dirty hack for loading iteration_results
    template<class T>
    void new_collector()
    {
        collector.reset(new detail::collector_impl<T>());
    }

    template<class T>
    void operator>>(T const& val)
    {
        if (!collector)
            collector.reset(new detail::collector_impl<T>());

    }

    const std::vector<boost::any>& get() const
    {
        return collector->get();
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
        for (std::map<std::string, boost::shared_ptr< ::detail::collector_impl_base> >::const_iterator
             it = collection.begin(); it != collection.end(); ++it)
        {
            ar[it->first] << *it->second;
        }
    }

    template <class Archive>
    void load(Archive & ar)
    {

        // TODO: This is dirty as hell
        // TODO: We must check the types of what comes out of the archive with the operator>>. Is there a way to do that?
        // For now, we check it manually
        std::vector<std::string> st = ar.list_children("");
        for (auto&& s: st)
        {
            // Create collectors beforehand
            if (s == "BondDimension") // for BondDimension use std::size_t
            {
                (*this)[s].new_collector<std::size_t>();
            }
            else // double
            {
                (*this)[s].new_collector<double>();
            }

            ar[s] >> *(collection[s].get());
        }


    }

    bool empty() const { return collection.empty(); };

private:
    std::map<std::string, boost::shared_ptr< ::detail::collector_impl_base> > collection;
};

#endif
