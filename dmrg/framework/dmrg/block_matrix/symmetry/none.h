/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#ifndef SYMMETRY_NONE_H
#define SYMMETRY_NONE_H

#include <iostream>

#include <boost/functional/hash.hpp>
 
#include <alps/hdf5.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/array.hpp>

class TrivialGroup
{
public:
	typedef enum { Plus } charge;
    typedef int subcharge; // Used if charge is site_dependent
	static const charge IdentityCharge = Plus;
    static const bool finite = true;
    
    
	static inline charge fuse(charge a, charge b) { return Plus; }
	template<int R> static charge fuse(boost::array<charge, R>) { return Plus; }
};

namespace boost {
    template <>
    class hash<TrivialGroup::charge>{
        public :
            size_t operator()(TrivialGroup::charge const &Charge ) const {
                return hash<int>()(TrivialGroup::IdentityCharge); // return 0 ??
            }
    };

    template <>
    class hash<std::pair<TrivialGroup::charge,TrivialGroup::charge > >{
        public :
            size_t operator()(std::pair<TrivialGroup::charge, TrivialGroup::charge> const &Pair_of_charge ) const {
                return boost::hash_value(Pair_of_charge);
            }
    };
};

inline void save(alps::hdf5::archive & ar
                 , std::string const & path
                 , TrivialGroup::charge const & value
                 , std::vector<std::size_t> size = std::vector<std::size_t>()
                 , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                 , std::vector<std::size_t> offset = std::vector<std::size_t>())
{
    int k = 1;
    ar[path] << k;
}

inline void load(alps::hdf5::archive & ar
                 , std::string const & path
                 , TrivialGroup::charge & value
                 , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                 , std::vector<std::size_t> offset = std::vector<std::size_t>())
{             
    value = TrivialGroup::Plus;
}

template <class Archive>
inline void serialize(Archive & ar, TrivialGroup::charge & c, const unsigned int version)
{
    ar & c;
}

inline TrivialGroup::charge operator-(TrivialGroup::charge a) { return a; }

#endif
