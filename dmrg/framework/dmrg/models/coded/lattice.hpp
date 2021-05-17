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

#ifndef MAQUIS_DMRG_lattice_hpp
#define MAQUIS_DMRG_lattice_hpp

#include "dmrg/models/lattice.h"

#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include<numeric>

#include "dmrg/utils/BaseParameters.h"

class ChainLattice : public lattice_impl
{
public:
    typedef lattice_impl::pos_t pos_t;

    ChainLattice (BaseParameters & parms, bool pbc_=false)
    : L(parms["L"])
    , a(parms["a"])
    , pbc(pbc_)
    { }

    ChainLattice (int L_, bool pbc_=false, double a_=1.)
    : L(L_)
    , a(a_)
    , pbc(pbc_)
    { }

    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (pbc && i == L-1)
            ret.push_back(0);
        return ret;
    }
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        if (pbc && i == L-1)
            ret.push_back(0);
        if (pbc && i == 0)
            ret.push_back(L-1);
        return ret;
    }

    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( site_label(pos[0]) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(pos[0], pos[1]) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( 0 );
        else if (property == "type" && pos.size() == 2)
            return boost::any( 0 );
        else if (property == "x" && pos.size() == 1)
            return boost::any( a * pos[0] );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L-1 );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( (pos[0] < pos[1]) );
        else if (property == "NumTypes")
            return boost::any( 1 );
        else if (property == "ParticleType" && pos.size() == 1)
            return boost::any( 0 );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    pos_t size() const
    {
        return L;
    }

    int maximum_vertex_type() const
    {
        return 0;
    }

private:

    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(a * i) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(a * i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(a * j) + " )");
    }

private:
    int L;
    double a;
    bool pbc;

};

class Orbitals : public lattice_impl
{
    typedef int subcharge;

public:
    typedef lattice_impl::pos_t pos_t;

    Orbitals (BaseParameters & parms)
    : L(parms["L"])
    , irreps(L, 0)
    , order(L)
    {
        if (!parms.is_set("orbital_order"))
            for (pos_t p = 0; p < L; ++p)
                order[p] = p;
        else {
            order = parms["orbital_order"].as<std::vector<pos_t> >();
            if (order.size() != L)
                throw std::runtime_error("Number of orbitals in the orbital order does not match the total number of orbitals");

            for (auto&& o: order) o--;
        }

        if (parms.is_set("integral_file")) {
            std::string integral_file = parms["integral_file"];
            if (!boost::filesystem::exists(integral_file))
                throw std::runtime_error("integral_file " + integral_file + " does not exist\n");

            std::ifstream orb_file;
            orb_file.open(parms["integral_file"].c_str());
            orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

            std::string line;
            std::getline(orb_file, line);

            orb_file.close();

            std::vector<std::string> split_line;
            boost::split(split_line, line, boost::is_any_of("="));

            // record the site_types in parameters
            parms.set("site_types", split_line[1]);
            irreps = parse_irreps(split_line[1]);
        }
        else if (parms.is_set("site_types"))
        {
            std::vector<subcharge> symm_vec = parms["site_types"];

            assert(L == symm_vec.size());
            for (subcharge p = 0; p < L; ++p)
                irreps[p] = symm_vec[order[p]];
        }
        else
            throw std::runtime_error("\"integral_file\" in parms input file or site_types is not set\n");

        maximum_vertex = *std::max_element(irreps.begin(), irreps.end());
    }

    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);

        return ret;
    }

    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        return ret;
    }

    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1) // return "( label )" as string
            return boost::any( site_label(order[pos[0]]) );
        if (property == "label_int" && pos.size() == 1) // just return label as integer
            return boost::any(order[pos[0]]);
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(order[pos[0]], order[pos[1]]) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( irreps[pos[0]] );
        else if (property == "type" && pos.size() == 2)
            return boost::any( 0 );
        else if (property == "NumTypes")
            return boost::any( 1 );
        else if (property == "ParticleType" && pos.size() == 1)
            return boost::any( 0 );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    pos_t size() const
    {
        return L;
    }

    int maximum_vertex_type() const
    {
        return maximum_vertex;
    }

private:

    std::vector<subcharge> parse_irreps(std::string input)
    {
        std::vector<subcharge> symm_vec, ret(L, 0);

        std::replace(input.begin(), input.end(), ',', ' ');
        std::istringstream iss(input);
        subcharge number;
        while( iss >> number )
            symm_vec.push_back(number-1);

        assert(L == symm_vec.size());
        for (subcharge p = 0; p < L; ++p)
            ret[p] = symm_vec[order[p]];

        return ret;
    }

    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(j) + " )");
    }

private:
    pos_t L;
    int maximum_vertex;

    std::vector<subcharge> irreps;
    std::vector<pos_t> order;
};

class SquareLattice : public lattice_impl
{
public:
    SquareLattice(BaseParameters & parms)
    : L_(parms["L"])
    , W_(parms["W"])
    , a(parms["a"])
    { }

    /*
     0 4  8 12
     1 5  9 13
     2 6 10 14
     3 7 11 15
     */
    std::vector<int> forward(int p) const
    {
        std::vector<int> ret;
        if (p+1 < L_*W_ && (p+1) % W_ != 0)
            ret.push_back(p+1);
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);

        //        maquis::cout << p << " -> ";
        //        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(maquis::cout, " "));
        //        maquis::cout << std::endl;

        return ret;
    }

    std::vector<int> all(int p) const
    {
        std::vector<int> ret = forward(p);
        if (p >= 1 && p % W_ != 0)
            ret.push_back(p-1);
        if (p >= W_)
            ret.push_back(p-W_);

        return ret;
    }

    int size() const { return L_*W_; }

    int maximum_vertex_type() const
    {
        return 0;
    }


    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( site_label(pos[0]) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(pos[0], pos[1]) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( 0 );
        else if (property == "type" && pos.size() == 2)
            return boost::any( 0 );
        else if (property == "x" && pos.size() == 1)
            return boost::any( x(pos[0]) );
        else if (property == "y" && pos.size() == 1)
            return boost::any( y(pos[0]) );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( false );
        else if (property == "NumTypes")
            return boost::any( 1 );
        else if (property == "ParticleType")
            return boost::any( 0 );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

private:

    double x (int i) const
    { return a * int(i/W_); }
    double y (int i) const
    { return a * (i%W_); }

    std::string site_label (int i) const
    {
        return "( " + ( boost::lexical_cast<std::string>(x(i))
                       + "," + boost::lexical_cast<std::string>(y(i)) ) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + ( boost::lexical_cast<std::string>(x(i))
                          + "," + boost::lexical_cast<std::string>(y(i)) ) + " )"
                + " -- "
                + "( " + ( boost::lexical_cast<std::string>(x(j))
                          + "," + boost::lexical_cast<std::string>(y(j)) ) + " )" );
    }


    int L_, W_;
    double a;
};


// +----------------------------+
//  PRE BORN OPPENHEIMER LATTICE
// +----------------------------+
//
// The PreBOLattice describes the sites for the full Hamiltonian.
//
class PreBOLattice : public lattice_impl
{
public:
    // Types definition
    typedef lattice_impl::pos_t  pos_t;
    // -- Constructor --
    // In addition to a standard lattice constructor, it also loads the number
    // of basis function per mode
    explicit PreBOLattice (BaseParameters & parms)
    {

        // Populate variables:
        auto vec_particles_str = parms["PreBO_ParticleTypeVector"].as<std::string>();
        auto isFermion_str = parms["PreBO_FermionOrBosonVector"].as<std::string>();
        auto orbitals_str = parms["PreBO_OrbitalVector"].as<std::string>();
        auto vec_ini_state_str = parms["PreBO_InitialStateVector"].as<std::string>();
        std::string max_m_str;
        if (parms.is_set("PreBO_MaxBondDimVector"))
            max_m_str = parms["PreBO_MaxBondDimVector"].as<std::string>();
        // convert strings to vectors
        std::istringstream is( vec_particles_str );
        vec_particles.assign(std::istream_iterator<int>( is ), std::istream_iterator<int>() );
        is.str(std::string());
        is.clear();
        is.str(isFermion_str);
        isFermion.assign(std::istream_iterator<int>( is ), std::istream_iterator<int>() );
        is.str(std::string());
        is.clear();
        is.str(orbitals_str);
        vec_orbitals.assign(std::istream_iterator<int>( is ), std::istream_iterator<int>() );
        is.str(std::string());
        is.clear();
        is.str(vec_ini_state_str);
        vec_ini_state.assign(std::istream_iterator<int>( is ), std::istream_iterator<int>() );
        is.str(std::string());
        is.clear();
        if (parms.is_set("PreBO_MaxBondDimVector")) {
            is.str(max_m_str);
            vec_max_m.assign(std::istream_iterator<int>( is ), std::istream_iterator<int>() );
            is.str(std::string());
            is.clear();
        }

        num_particle_types = vec_particles.size();
        // ATTENTION MINUS ONE!!!
        maximum_vertex = num_particle_types-1;

        // construct lattice containing the particle types
        // the length of the lattice is determined by the number of different particle types and the number
        // of the basis functions of each particle type.
        L = std::accumulate(vec_orbitals.begin(), vec_orbitals.end(), 0);
        parms.set("L", L);
        std::vector<part_type> vec_abs_index_part_type;
        for (part_type i = 0; i < num_particle_types; i++) {
            for (std::size_t j = 0; j < vec_orbitals[i]; j++) {
                vec_abs_index_part_type.push_back(i);
            }
        }
        if (parms.is_set("orbital_order")) {
            m_order = parms["orbital_order"].as<std::vector<pos_t> >();
            vec_lattice_type.resize(L);
            m_inv_order.resize(L);
            if (m_order.size() != L)
                throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");
            for (int p = 0; p < m_order.size(); ++p)
                m_inv_order[p] = std::distance(m_order.begin(), std::find(m_order.begin(), m_order.end(), p));
            for (pos_t i = 0; i < L; i++) {
                vec_lattice_type[i] = vec_abs_index_part_type[m_order[i]];
            }
        }
        else
            vec_lattice_type=vec_abs_index_part_type;
        //populate vec_fer_bos
        unsigned int bos_temp = 0;
        unsigned int fer_temp = 0;
        for (unsigned int i=0; i<num_particle_types; i++) {
            if (isFermion[i]) {
                vec_fer_bos.push_back(fer_temp);
                fer_temp++;
            }
            else {
                vec_fer_bos.push_back(bos_temp);
                bos_temp++;
            }
        } // for loop


    } // Constructor
    // +-------+
    //  Methods
    // +-------+
    // -- GETTERS --
    /**
     * Takes the relative position of a given particle type and returns its absolute position on the lattice
     * @param pt
     * @param rel_pos
     * @return
     */
    pos_t get_abs_position(part_type const & pt, pos_t const & rel_pos) const {
        //throw std::runtime_error("get_abs_position must be debugged first.");
        //unsigned int abs_pos=0;
        //for (unsigned int i=0; i<pt; i++) {
        //    for (unsigned int j=0; j<vec_orbitals[i]; j++) {
        //        abs_pos++;
        //    }
        //}
        //abs_pos+=rel_pos;
        //return abs_pos;
//        auto it = (vec_orbitals.begin() + pt);
        return std::accumulate(vec_orbitals.begin(), (vec_orbitals.begin() + pt) , rel_pos);
    }
    // The following methods are the same as for the Orbital and the Open Chain Lattice.
    // The only difference is the way in which the lattice site type is extracted, which
    // is taken from the orbital lattice
    std::vector<pos_t> forward(pos_t i) const {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        return ret;
    }
    //
    std::vector<pos_t> all(pos_t i) const {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        return ret;
    }
    //
    /**
     * This function takes the property as input and returns its value.
     * E.g.: "type" at position 12 in lattice --> 0 (electron)
     * @param property
     * @param pos
     * @return boost::any
     */
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "type" && pos.size() == 1)
            return boost::any(vec_lattice_type[pos[0]]);
        else if (property == "Mmax" && pos.size() == 1)
            return boost::any(vec_max_m[pos[0]]);
        else if (property == "NumTypes")
            return boost::any(num_particle_types);
        else if (property == "ParticleType" && pos.size() == 1)
            return boost::any( vec_lattice_type[pos[0]] );
        else if (property == "label" && pos.size() == 2)
            return boost::any(bond_label(pos[0], pos[1]));
        else if (property == "vec_particles")
            return boost::any(vec_particles);
        else if (property == "num_particle_types")
            return boost::any(num_particle_types);
        else if (property == "isFermion")
            return boost::any(isFermion);
        else if (property == "vec_orbitals")
            return boost::any(vec_orbitals);
        else if (property == "vec_ini_state")
            return boost::any(vec_ini_state);
        else if (property == "vec_fer_bos")
            return boost::any(vec_fer_bos);
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
            throw std::runtime_error(ss.str());
        }
    }

    pos_t size()                const { return L; }
    int   maximum_vertex_type() const { return maximum_vertex; }

private:
    // +----------+
    //  Attributes
    // +----------+
    pos_t L;       // Size of the DMRG lattice
    int num_particle_types;                    // Number of particle types
    std::vector<int> vec_particles;            // Number of particles per type
    std::vector<int> vec_ini_state;            // Specifies the initial state of the system.
    std::vector<bool> isFermion;               // Fermion 1, Boson 0
    std::vector<int> vec_orbitals;             // Number of orbitals per type
    std::vector<part_type> vec_lattice_type;   // Stores the index of the particle type for each site.
    int maximum_vertex;                        // Largest index for site types
    std::vector<int> vec_ibs;                  // Just in here for safety checks
    std::vector<int> vec_fer_bos;              // vector that maps the particle types vector
    std::vector<pos_t> m_order;                // ordering of the sites -> E.g. canonical or Fiedler
    std::vector<pos_t> m_inv_order;            // inverted ordering of the sites
    std::vector<std::size_t> vec_max_m;        // Stores max bond dim for all particle types

    // +-------------------+
    //   Printing routines
    // +-------------------+
    std::string site_label (int i) const {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }

    std::string bond_label (int i, int j) const {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
                  + " -- " + "( " + boost::lexical_cast<std::string>(j) + " )");
    }
} ;

#endif
