/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Michele Dolfi <dolfim@phys.ethz.ch>
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
