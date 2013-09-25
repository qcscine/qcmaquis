/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef CI_DEAS_HPP
#define CI_DEAS_HPP

#include <iostream>

#include <vector>
#include <set>
#include <map>
#include <string>
#include <boost/lexical_cast.hpp>

#include "dmrg/utils/DmrgParameters2.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_initializers.h"


std::vector<std::string>
parse_config_strings(std::string file)
{
    std::ifstream config_file;
    config_file.open(file.c_str());

    std::vector<std::string> configs;

    for (std::string line; std::getline(config_file, line); ) {
        std::vector<std::string> det_raw;
        boost::split(det_raw, line, boost::is_any_of(" "));
        
        std::string det = det_raw[0];
        maquis::cout << "parsed " << det << std::endl;
        configs.push_back(det);
    }
    
    return configs;
}

std::vector< std::map<typename TwoU1::charge, std::set<std::string> > >
arrange_configs(std::vector<std::string> input_configs)
{
    typedef TwoU1 grp;
    typedef typename grp::charge charge;
    typedef unsigned pos_t;

    TwoU1::charge A(1), B(0), C(0), D(0);
    B[0]=1; C[1]=1;

    pos_t config_length = input_configs[0].size();
    maquis::cout << "config length " << config_length << std::endl;
    maquis::cout << "nconfigs " << input_configs.size() << std::endl;

    std::vector< std::map<charge, std::set<std::string> > > ret(config_length);

    for (std::size_t c = 0; c < input_configs.size(); ++c)
    {
        maquis::cout << "follow config " << c << std::endl;
        std::string const & config = input_configs[c];
        charge sector = grp::IdentityCharge;
        for (int p = config_length-1; p >= 0; --p)
        {
            int occ = boost::lexical_cast<int>(config[p]);
            charge onsite;
            switch(occ) {
                case 4:
                    onsite = A;
                    break;
                case 3:
                    onsite = B;
                    break;
                case 2:
                    onsite = C;
                    break;
                case 1:
                    onsite= D;
                    break;
            }

            sector = grp::fuse(sector, onsite);
            //std::pair< std::set<std::string>::const_iterator, bool > iinfo;
            ret[p][sector].insert(config.substr(p)); 
        }
    }
    return ret;
}

void display_environment(std::vector< std::map<typename TwoU1::charge, std::set<std::string> > > const & env)
{
    typedef TwoU1 grp;
    typedef typename grp::charge charge;

    for (std::size_t p = 0; p < env.size(); ++p)
    {
        maquis::cout << "env position " << p << std::endl << std::endl;
        for (typename std::map<charge, std::set<std::string> >::const_iterator itm = env[p].begin();
                itm != env[p].end(); ++itm)
        {
            maquis::cout << "charge sector " << itm->first << std::endl;
            for (typename std::set<std::string>::const_iterator its = itm->second.begin();
                    its != itm->second.end(); ++its) {
                maquis::cout << *its << std::endl;
            }
            maquis::cout << std::endl;
        }
        maquis::cout << std::endl;
    }
}

#endif
