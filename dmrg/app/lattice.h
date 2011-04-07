/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef LATTICE_H
#define LATTICE_H

#include <vector>
#include <string>
#include <fstream>
#include <boost/any.hpp>

class Lattice
{
public:
    typedef int pos_t;
    
    virtual std::vector<pos_t> forward(pos_t) const = 0;
    virtual std::vector<pos_t> all(pos_t) const = 0;
    
    // non-virtual!
    template<class T> T get_prop(std::string property,
                                 pos_t site) const
    {
        return boost::any_cast<T>(get_prop_(property, std::vector<pos_t>(1, site)));
    }
    
    template<class T> T get_prop(std::string property,
                                 pos_t bond1, pos_t bond2) const
    {
        std::vector<pos_t> v(2);
        v[0] = bond1; v[1] = bond2;
        return boost::any_cast<T>(get_prop_(property, v));
    }
    
    template<class T> T get_prop(std::string property,
                                 std::vector<pos_t> const & positions) const
    {
        return boost::any_cast<T>(get_prop_(property, positions));
    }
    
    // virtual!
    virtual boost::any get_prop_(std::string const &, std::vector<pos_t> const &) const = 0;
    
    virtual pos_t size() const = 0;
};


Lattice * lattice_factory (std::string const & lattice, std::ifstream & ifs);
#endif
