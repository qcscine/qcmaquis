//
//  cont_lattice.h
//  
//
//  Created by Michele Dolfi on 30.06.11.
//  Copyright 2011 ETHZ. All rights reserved.
//

#ifndef CONTINUOUS_LATTICE_H
#define CONTINUOUS_LATTICE_H

#include "lattice.h"

#include <sstream>
#include <boost/lexical_cast.hpp>

#include "utils/DmrgParameters.h"

namespace app {
    
    class ContChain : public Lattice
    {
    public:
        typedef Lattice::pos_t pos_t;
        
        ContChain (ModelParameters & parms, bool pbc_=false)
        : L(parms.get<int>("L"))
        , N(parms.get<int>("Ndiscr"))
        , a(parms.get<double>("a"))
        , pbc(pbc_)
        {}
        
        std::vector<pos_t> forward(pos_t i) const
        {
            std::vector<pos_t> ret;
            if (i < L*N-1)
                ret.push_back(i+1);
            if (pbc && i == L*N-1)
                ret.push_back(0);
            return ret;
        }
        std::vector<pos_t> all(pos_t i) const
        {
            std::vector<pos_t> ret;
            if (i < L*N-1)
                ret.push_back(i+1);
            if (i > 0)
                ret.push_back(i-1);
            if (pbc && i == L*N-1)
                ret.push_back(0);
            if (pbc && i == 0)
                ret.push_back(L*N-1);
            return ret;
        }
        
        boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
        {
            if (property == "label" && pos.size() == 1)
                return boost::any( site_label(pos[0]) );
            else if (property == "label" && pos.size() == 2)
                return boost::any( bond_label(pos[0], pos[1]) );
            else if (property == "type" && pos.size() == 1)
                return boost::any( pos[0]%N );
            else if (property == "type" && pos.size() == 2)
                return boost::any( 0 );
            else if (property == "x" && pos.size() == 1)
                return boost::any( a/N * pos[0] );
            else if (property == "dx" && pos.size() == 1)
                return boost::any( a/N );
            else if (property == "dx" && pos.size() == 2)
                return boost::any( a/N );
            else if (property == "wraps_pbc" && pos.size() == 2)
                return boost::any( (pos[0] < pos[1]) );
            else {
                std::ostringstream ss;
                ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
                throw std::runtime_error(ss.str());
                return boost::any();
            }
        }
        
        pos_t size() const
        {
            return L*N;
        }
        
    private:
        
        std::string site_label (int i) const
        {
            return "( " + boost::lexical_cast<std::string>(a/N * i) + " )";
        }
        
        std::string bond_label (int i, int j) const
        {
            return (  "( " + boost::lexical_cast<std::string>(a/N * i) + " )"
                    + " -- "
                    + "( " + boost::lexical_cast<std::string>(a/N * j) + " )");
        }
        
    private:
        int N, L;
        double a;
        bool pbc;
        
    };
    
} // namespace

#endif
