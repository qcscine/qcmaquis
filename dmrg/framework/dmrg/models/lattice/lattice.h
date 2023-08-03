/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef LATTICE_H
#define LATTICE_H

#include "dmrg/utils/BaseParameters.h"

#include <utility>
#include <vector>
#include <string>
#include <boost/any.hpp>

/// lattice common base
class lattice_impl {
public:
    using pos_t = int;
    using part_type = int;

    virtual ~lattice_impl() = default;

    virtual std::vector<pos_t> forward(pos_t) const = 0;
    virtual std::vector<pos_t> all(pos_t) const = 0;

    // non-virtual!
    template<class T> T get_prop(const std::string& property) const
                                 
    {
        return boost::any_cast<T>(get_prop_(property, std::vector<pos_t>()));
    }

    template<class T> T get_prop(const std::string& property,
                                 pos_t site) const
    {
        return boost::any_cast<T>(get_prop_(property, std::vector<pos_t>(1, site)));
    }

    template<class T> T get_prop(const std::string& property,
                                 pos_t bond1, pos_t bond2) const
    {
        std::vector<pos_t> v(2);
        v[0] = bond1; v[1] = bond2;
        return boost::any_cast<T>(get_prop_(property, v));
    }

    template<class T> T get_prop(const std::string& property,
                                 std::vector<pos_t> const & positions) const
    {
        return boost::any_cast<T>(get_prop_(property, positions));
    }

    // virtual!
    virtual boost::any get_prop_(std::string const &, std::vector<pos_t> const &) const = 0;

    virtual pos_t get_abs_position(part_type const & pt, pos_t const & rel_pos) const {return 0;};

    virtual pos_t size() const = 0;
    
    /** @brief Getter for the number of types available in the lattice */
    virtual int getMaxType() const = 0;

};


/// lattice factory
std::shared_ptr<lattice_impl>
lattice_factory(BaseParameters & parms);


/// pimpl resolved Lattice
class Lattice {
    using impl_type = lattice_impl;
    using impl_ptr = std::shared_ptr<lattice_impl>;
public:
    using pos_t = impl_type::pos_t;
    using part_type = int;

    Lattice() = default;

    Lattice(BaseParameters & parms)
    : impl_(lattice_factory(parms))
    { }

    Lattice(impl_ptr impl) : impl_(std::move(impl)) { }

    impl_ptr impl() const { return impl_; }

    std::vector<pos_t> forward(pos_t site) const { return impl_->forward(site); }
    std::vector<pos_t> all(pos_t site) const { return impl_->all(site); }
    
    /** @brief General lattice properties */
    template<class T> T get_prop(std::string property) const
    { return impl_->get_prop<T>(property); }

    /** @brief Site-specific lattice properties */
    template<class T> T get_prop(std::string property, pos_t site) const
    { return impl_->get_prop<T>(property, site); }

    /** @brief Bond-specific lattice properties */
    template<class T> T get_prop(std::string property, pos_t bond1, pos_t bond2) const
    { return impl_->get_prop<T>(property, bond1, bond2); }

    /** @brief Multi-site properties */
    template<class T> T get_prop(std::string property, std::vector<pos_t> const & positions) const
    { return impl_->get_prop<T>(property, positions); }

    pos_t get_abs_position(part_type const & pt, pos_t const & rel_pos) const
    { return impl_->get_abs_position(pt, rel_pos) ; }

    pos_t size() const { return impl_->size(); }

    /** @brief Getter for the number of types available in the lattice */
    int getMaxType() const { return impl_->getMaxType(); };

private:
    impl_ptr impl_;
};



#endif
