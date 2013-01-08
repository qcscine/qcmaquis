/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "runsim.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include <boost/function.hpp>

#include "symm_factory.h"

namespace maquis { namespace dmrg {
    
    void symm_factory(DmrgParameters & parms, ModelParameters & model)
    {
        std::map<std::string, boost::function<void (DmrgParameters & p, ModelParameters & m)> > factory_map;
        
        maquis::cout << "This binary contains symmetries: ";
#ifdef HAVE_TrivialGroup
        factory_map["none"] = run_sim<TrivialGroup>;
        maquis::cout << "none ";
#endif
#ifdef HAVE_U1
        factory_map["u1"] = run_sim<U1>;
        maquis::cout << "u1 ";
#endif
#ifdef HAVE_TwoU1
        factory_map["2u1"] = run_sim<TwoU1>;
        maquis::cout << "2u1 ";
#endif
#ifdef HAVE_Ztwo
        factory_map["Z2"] = run_sim<Ztwo>;
        maquis::cout << "Z2 ";
#endif
        maquis::cout << std::endl;
        
        std::string symm_name = parms.get<std::string>("symmetry");
        
        if (factory_map.find(symm_name) != factory_map.end())
            factory_map[symm_name](parms, model);
        else
            throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");
    }

} }
