/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <robinfe@phys.chem.ethz.ch>
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

#ifndef MAQUIS_DMRG_WATSON_LATTICE
#define MAQUIS_DMRG_WATSON_LATTICE

#ifdef DMRG_VIBRATIONAL

#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/lattice/lattice.h"

/**
 * @brief Watson lattice class.
 * 
 * Lattice representing the canonical quantization-based vibrational Hamiltonian.
 * Unlike the NU1 case, each site is mapped to a mode.
 */

class WatsonLattice : public lattice_impl
{
public:
    // Types definition
    using post_t = lattice_impl::pos_t;
    
    /**
     * @brief Class constructor for the lattice
     * @param parameters parameter container
     */
    WatsonLattice(BaseParameters& parameters) : L(parameters["L"]), maximum_vertex(1)
    { }
    
    /** @brief Returns the next position in the lattice */
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        return ret;
    }
    
    /** @brief Returns the neighbors of a given site */
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        return ret;
    }

    /**
     * @brief Getter for the property
     * 
     * Note that, in additional to the usual properties of a lattice, we code
     * the additional property "sublatticePos" which states where the sublattice
     * associated with a given mode is starting
     * 
     * @param property string identifier for the property
     * @param pos vector of positions 
     * @return boost::any requested property
     */
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any(site_label(pos[0]));
        else if (property == "label" && pos.size() == 2)
            return boost::any(bond_label(pos[0], pos[1]));
        else if (property == "type" && pos.size() == 1)
            return boost::any(0);
        else if (property == "type" && pos.size() == 2)
            return boost::any(0);
        else if (property == "ParticleType" && pos.size() == 1) {
            assert (pos[0] >= 0 && pos[0] < L);
            return 0;
        }
        else if (property == "NumTypes")
            return 1;
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }
    
    /** @brief Getter for the lattice size */
    pos_t size() const { return L; }

    /** @brief Getter for the number of types of sites */
    int maximum_vertex_type() const { return maximum_vertex; }

private:
    /** Size of the DMRG lattice */
    pos_t L;
    /** Largest index for site types */
    int maximum_vertex;
    
    /** @brief Prints the label of a given site */
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }
    
    /** @brief Prints the label of a given bond */
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(j) + " )");
    }
};

#endif // DMRG_VIBRATIONAL

#endif // MAQUIS_DMRG_NMODE_LATTICE