/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef UTILS_GUESS_SYMMETRY_H
#define UTILS_GUESS_SYMMETRY_H

#include "dmrg/utils/DmrgParameters.h"

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <map>
#include <string>

inline std::string guess_alps_symmetry(BaseParameters & parms)
{
    std::map<int, std::string> symm_names;
    symm_names[0] = "none";
    symm_names[1] = "u1";
    symm_names[2] = "2u1";

    int n=0;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
        boost::char_separator<char> sep(" ,");
        std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"].str();
        tokenizer qn_tokens(qn_string, sep);
        for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
            if (parms.defined(*it + "_total"))
                n += 1;
        }
    }
    return symm_names[n];
}

#endif
