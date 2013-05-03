/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef SYMMETRY_NONE_H
#define SYMMETRY_NONE_H

#include <iostream>

#include <boost/functional/hash.hpp>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#include <boost/serialization/serialization.hpp>

class TrivialGroup
{
public:
	typedef enum { Plus } charge;
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

#ifdef HAVE_ALPS_HDF5
namespace alps {
    namespace hdf5 {
        
        inline void save(
                         alps::hdf5::archive & ar
                         , std::string const & path
                         , TrivialGroup::charge const & value
                         , std::vector<std::size_t> size = std::vector<std::size_t>()
                         , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                         , std::vector<std::size_t> offset = std::vector<std::size_t>())
        {
            int k = 1;
            ar << alps::make_pvp(path, k);            
        }
        inline void load(
                         alps::hdf5::archive & ar
                         , std::string const & path
                         , TrivialGroup::charge & value
                         , std::vector<std::size_t> chunk = std::vector<std::size_t>()
                         , std::vector<std::size_t> offset = std::vector<std::size_t>())
        {
            value = TrivialGroup::Plus;
        }
        
    }
}
#endif

template <class Archive>
inline void serialize(Archive & ar, TrivialGroup::charge & c, const unsigned int version)
{
    ar & c;
}

inline TrivialGroup::charge operator-(TrivialGroup::charge a) { return a; }

#endif
