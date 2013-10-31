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
#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>

#include "symm_factory.h"

std::string guess_alps_symmetry(ModelParameters& model)
{
    std::map<int, std::string> symm_names;
    symm_names[0] = "none";
    symm_names[1] = "u1";
    symm_names[2] = "2u1";
    
    int n=0;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    if (model.defined("CONSERVED_QUANTUMNUMBERS")) {
        boost::char_separator<char> sep(" ,");
        std::string qn_string = model["CONSERVED_QUANTUMNUMBERS"];
        tokenizer qn_tokens(qn_string, sep);
        for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
            if (model.defined(*it + "_total"))
                n += 1;
        }
    }
    return symm_names[n];
}

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
        
        std::string symm_name;
        if (parms["model_library"] == "alps")
            symm_name = guess_alps_symmetry(model);
        else
            symm_name = parms["symmetry"].str();
        
        if (factory_map.find(symm_name) != factory_map.end())
            factory_map[symm_name](parms, model);
        else
            throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");
#ifdef AMBIENT
        ambient::sync();
#endif
    }

} }
