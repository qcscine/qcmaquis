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

#include <boost/any.hpp>
#include <vector>
#include <memory>

/**
 * @brief Class used to store the results of a generic sweep-based algorithm 
 * The object is basically a wrapper around a string --> vector<boost::any> object.
 * The string identifies the specific result of the simulation.
 * The vector<boost::any> object is represented, in practice, by a pointer to a 
 * [collector_impl_base] object. The latter is, in turn, a interface class
 * that in implemented, by the [collector_impl] class in the "<<" operator.
 */
class results_collector
{
private:
    /** @brief Virtual class representing an object that can hold generic values (see cpp for implementation) */
    class collector_impl_base;

    /**
     * @brief Actual implementation based on a boost::any vector 
     * Note that the template parameter T represents the type used for the casting of boost::any.
     */
    template <class T>
    class collector_impl;
public:

    /**
     * @brief Proxy class for a pointer to a [collector_impl_base] class.
     * Note that [collector_impl_base] is a purely virtual class, so an interface.
     * So this proxy class does not assume anything about its implementation.
     */
    class collector_proxy {
    typedef std::shared_ptr<results_collector::collector_impl_base> coll_type;
    public:

        /** @brief Class constructor, just stores the reerence to the pointer */
        collector_proxy(coll_type & collector_) : collector(collector_) { }

        /** @brief Loading method */
        template<class T>
        void operator<<(T const& val);

        // Needed for a dirty hack for loading iteration_results
        template<class T>
        void new_collector();

        /** @brief Reading method */
        template<class T>
        void operator>>(T const& val);

        /** @brief Gets the underlying boost::get object */
        const std::vector<boost::any>& get() const;

    private:
        coll_type& collector;
    };

    /** @brief Getter of a given result */
    collector_proxy operator[] (std::string name);

    /** @brief Reset method */
    void clear();

    /** @brief Stores data in a hdf5 archive */
    template <class Archive>
    void save(Archive & ar) const;

    /** @brief Loads data from a hdf5 archive */
    template <class Archive>
    void load(Archive & ar);

    /** @brief Checks if a collector is empty */
    bool empty() const;

    /** @brief Checks the presence of a given element in the collector */
    bool has(const std::string key) const { return collection.find(key) != collection.end(); }

private:
    std::map<std::string, std::shared_ptr<collector_impl_base> > collection;
};

#endif
