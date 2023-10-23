/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include <map>
#include "results_collector.h"
#include "dmrg/utils/storage.h"
#include <boost/preprocessor/seq/for_each.hpp>
#include "dmrg/sim/matrix_types.h"

class results_collector::collector_impl_base
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
class results_collector::collector_impl : public results_collector::collector_impl_base
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

// results_collector::collector_proxy implementation

template<class T>
void results_collector::collector_proxy::operator<<(T const& val)
{
    if (!collector)
        collector.reset(new results_collector::collector_impl<T>());
    collector->collect(val);
}

template<class T>
void results_collector::collector_proxy::new_collector()
{
    collector.reset(new results_collector::collector_impl<T>());
}

template<class T>
void results_collector::collector_proxy::operator>>(T const& val)
{
    if (!collector)
        collector.reset(new results_collector::collector_impl<T>());
}

const std::vector<boost::any>& results_collector::collector_proxy::get() const
{
    return collector->get();
}

// ------ results_collector implementation --------
void results_collector::clear()
{
    collection.clear();
}

results_collector::collector_proxy results_collector::operator[] (std::string name)
{
    return results_collector::collector_proxy(collection[name]);
}

template <class Archive>
void results_collector::save(Archive & ar) const
{
    for (const auto& it : collection)
    {
        ar[it.first] << *it.second;
    }
}

template <class Archive>
void results_collector::load(Archive & ar)
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

bool results_collector::empty() const { return collection.empty(); };

// instantiate template functions
template void results_collector::save<alps::hdf5::archive>(alps::hdf5::archive&) const;
template void results_collector::load<alps::hdf5::archive>(alps::hdf5::archive&);

#define INSTANTIATE_COLLECTOR_PROXY(r, d, T) \
template void results_collector::collector_proxy::operator<< <T>(T const&);

#define COLLECTOR_PROXY_TYPES (matrix::value_type) (cmatrix::value_type) (unsigned long)

BOOST_PP_SEQ_FOR_EACH(INSTANTIATE_COLLECTOR_PROXY, _, COLLECTOR_PROXY_TYPES)

#undef INSTANTIATE_COLLECTOR_PROXY
#undef COLLECTOR_PROXY_TYPES