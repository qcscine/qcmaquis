/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

#ifndef MAQUIS_DMRG_NMODE_LATTICE
#define MAQUIS_DMRG_NMODE_LATTICE

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"

// +-------------+
//  NMODE LATTICE
// +-------------+
//TODO ALB Maybe move the setup of the orbital order here. When doing it, remember to change also get_prop_

class NModeLattice : public lattice_impl
{
public:
    // Types definition
    typedef lattice_impl::pos_t   pos_t;
    //
    // -- Constructor --
    // In addition to a standard lattice constructor, it also loads the number
    // of basis function per mode
    NModeLattice (BaseParameters & model) : L(0), vector_bases(0), vector_types(0)
    {
        // Loads in input parameters
        int num_modes = model["nmode_num_modes"] ;
        L = model["L"] ;
        vector_types.resize(L) ;
        vector_bases.reserve(L) ;
        std::string jnk = model["nmode_num_basis"];
        maximum_vertex = num_modes-1 ;
        std::vector<std::string>  size_vec ;
        boost::split(size_vec, jnk, boost::is_any_of(",")) ;
        for (std::size_t idx = 0; idx < size_vec.size(); idx++) {
            if (idx == 0) {
                vector_bases.push_back(0);
            } else {
                int mod = stoi(size_vec[idx-1]);
                if (mod <= 0)
                    throw std::runtime_error("Non-positive number of basis function found");
                else
                    vector_bases.push_back(vector_bases[idx-1] + mod);
            }
        }
        int count = vector_bases[size_vec.size()-1] + stoi(size_vec[size_vec.size()-1])  ;
        if (count != L)
            throw std::runtime_error("Inconsistent number of basis functions") ;
        // Now populates the vector with the type of each site (i.e., the mode to
        // which they belong
        std::size_t jcont=0 ;
        for (std::size_t idx1 = 0; idx1 < size_vec.size(); idx1++) {
            for (std::size_t idx2 = 0; idx2 < stoi(size_vec[idx1]); idx2++) {
                vector_types[jcont] = static_cast<int>(idx1) ;
                ++jcont ;
            }
        }
    }
    //
    // -- METHODS --
    // The following methods are the same as for the Orbital and the Open Chain Lattice.
    // The only difference is the way in which the lattice site type is extracted, which
    // is taken from the orbital lattice
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        return ret;
    }
    //
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        return ret;
    }
    //
    int get_basis_previous(int const& idx) const
    {
        assert (idx >= 0 && idx < L) ;
        return vector_bases[idx] ;
    }

    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any(site_label(pos[0]));
        else if (property == "label" && pos.size() == 2)
            return boost::any(bond_label(pos[0], pos[1]));
        else if (property == "type" && pos.size() == 1)
            return boost::any(vector_types[pos[0]]);
        else if (property == "type" && pos.size() == 2)
            return boost::any(0);
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }
    //
    pos_t size()                const { return L; } 
    int   maximum_vertex_type() const { return maximum_vertex; }
private:
    //
    // -- ATTRIBUTES --
    //
    pos_t L;                                  // Size of the DMRG lattice
    int maximum_vertex;                       // Largest index for site types
    std::vector<int> vector_types ;           // Sites type vector
    std::vector<int> vector_bases ;           // The i-th element returns the number of basis that have been used
                                              // before the i-th mode. Used as offset in vectors
    //
    // Printing routines
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }
    //
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(j) + " )");
    }
};

#endif // MAQUIS_DMRG_NMODE_LATTICE