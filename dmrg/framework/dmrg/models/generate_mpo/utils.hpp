/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef GENERATE_MPO_UTILS_H
#define GENERATE_MPO_UTILS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/models/OperatorHandlers/OpTable.h"
#include "dmrg/models/lattice/lattice.h"

#include <string>
#include <sstream>

#include <boost/bind.hpp>

namespace generate_mpo
{
	template<class Matrix, class SymmGroup>
	struct OperatorTagTerm
	{
		typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Lattice::pos_t pos_t;
		typedef std::pair<pos_t, tag_type> op_pair_t;
        
		std::vector<op_pair_t> operators;
		tag_type fill_operator;
        typename Matrix::value_type scale;
        bool with_sign;
        
        OperatorTagTerm() : scale(1.), with_sign(false) {}
        
        void canonical_order() // TODO: check and fix for fermions
        {
            std::sort(operators.begin(), operators.end(),
                      boost::bind(&op_pair_t::first, _1) <
                      boost::bind(&op_pair_t::first, _2));
        }
        
        bool operator< (OperatorTagTerm const & rhs) const
        {
            if (operators[0].first == rhs.operators[0].first)
                return operators.size() >= rhs.operators.size();
            return operators[0].first < rhs.operators[0].first;
        }
        
        bool site_match (OperatorTagTerm const & rhs) const
        {
            if (operators.size() == rhs.operators.size())
            {
                bool ret = true;
                for (std::size_t p=0; p<operators.size() && ret; ++p)
                    ret = (operators[p].first == rhs.operators[p].first);
                return ret;
            } else if (operators.size() == 2 && rhs.operators.size() == 1)
                return (operators[0].first == rhs.operators[0].first || operators[1].first == rhs.operators[0].first);
            else if (operators.size() == 1 && rhs.operators.size() == 2)
                return (operators[0].first == rhs.operators[0].first || operators[0].first == rhs.operators[1].first);
            else
            {
                throw std::runtime_error("site_match not implemented for this type of operator." );
                return false;
            }
            
        }
        
        bool overlap (OperatorTagTerm const & rhs) const
        {
        	return !( (operators.rbegin()->first < rhs.operators.begin()->first) || (rhs.operators.rbegin()->first < operators.begin()->first) );
        }
	};
    
    template<class Matrix, class SymmGroup>
    std::ostream & operator<< (std::ostream & os, OperatorTagTerm<Matrix, SymmGroup> const& op)
    {
        os << "fill: " << op.fill_operator << std::endl;
        os << "sign: " << op.with_sign << std::endl;
        os << "scale: " << op.scale << std::endl;
        os << "operators:";
        for (int i=0; i<op.operators.size(); ++i)
            os << " {"  << op.operators[i].first << "," << op.operators[i].second << "}";
            os << std::endl;
        return os;
    }
    
	template<class Matrix, class SymmGroup>
	struct OperatorTerm
	{
		typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef Lattice::pos_t pos_t;
		typedef std::pair<pos_t, op_t> op_pair_t;
        
		std::vector<op_pair_t> operators;
		op_t fill_operator;
        bool with_sign;
        
        OperatorTerm() : with_sign(false) {}
        
        void canonical_order() // TODO: check and fix for fermions
        {
            std::sort(operators.begin(), operators.end(),
                      boost::bind(&op_pair_t::first, _1) <
                      boost::bind(&op_pair_t::first, _2));
        }
        
        bool operator< (OperatorTerm const & rhs) const
        {
            if (operators[0].first == rhs.operators[0].first)
                return operators.size() >= rhs.operators.size();
            return operators[0].first < rhs.operators[0].first;
        }

        bool site_match (OperatorTerm const & rhs) const
        {
            if (operators.size() == rhs.operators.size())
            {
                bool ret = true;
                for (std::size_t p=0; p<operators.size() && ret; ++p)
                    ret = (operators[p].first == rhs.operators[p].first);
                return ret;
            } else if (operators.size() == 2 && rhs.operators.size() == 1)
                return (operators[0].first == rhs.operators[0].first || operators[1].first == rhs.operators[0].first);
            else if (operators.size() == 1 && rhs.operators.size() == 2)
                return (operators[0].first == rhs.operators[0].first || operators[0].first == rhs.operators[1].first);
            else
            {
                throw std::runtime_error("site_match not implemented for this type of operator." );
                return false;
            }
                
        }
        
        bool overlap (OperatorTerm const & rhs) const
        {
        	return !( (operators.rbegin()->first < rhs.operators.begin()->first) || (rhs.operators.rbegin()->first < operators.begin()->first) );
        }

	};
   
    using namespace std;
    using namespace boost::tuples;

    inline size_t next_free(vector<size_t> const & out_taken,
                            vector<size_t> const & in_taken)
    {
        for (size_t k = 0; true; ++k)
        {
            if (count(out_taken.begin(), out_taken.end(), k) == 0 &&
                count(in_taken.begin(), in_taken.end(), k) == 0)
                return k;// +-------------+
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
        }
    }
    
    inline size_t next_free(set<size_t> const & s)
    {
        for (size_t k = 2; true; ++k)
            if (s.count(k) == 0)
                return k;
    }
    
    template<class Vector>
    void compress_on_bond(Vector & pm1, Vector & pm2)
    {
        std::set<size_t> bond_used_dims;
        for (typename Vector::iterator it = pm1.begin(); it != pm1.end(); ++it)
            if (get<1>(*it) > 1)
                bond_used_dims.insert(get<1>(*it));
        for (typename Vector::iterator it = pm2.begin(); it != pm2.end(); ++it)
            if (get<0>(*it) > 1)
                bond_used_dims.insert(get<0>(*it));
        
        std::map<size_t, size_t> compression_map;
        size_t c = 2;
        for (set<size_t>::iterator it = bond_used_dims.begin();
             it != bond_used_dims.end(); ++it)
            compression_map[*it] = c++;
        
        for (typename Vector::iterator it = pm1.begin(); it != pm1.end(); ++it)
            if (compression_map.count(get<1>(*it)) > 0)
                get<1>(*it) = compression_map[get<1>(*it)];
        for (typename Vector::iterator it = pm2.begin(); it != pm2.end(); ++it)
            if (compression_map.count(get<0>(*it)) > 0)
                get<0>(*it) = compression_map[get<0>(*it)];
    }

    template<class Vector>
    std::pair<size_t, size_t> rcdim(Vector const & pm)
    {
        std::list<size_t> l, r;
        for (typename Vector::const_iterator it = pm.begin(); it != pm.end(); ++it) {
            l.push_back( get<0>(*it) );
            r.push_back( get<1>(*it) );
        }
        
        size_t ldim=0, rdim=0;
        if (l.size() > 0) ldim = *max_element(l.begin(), l.end())+1;
        if (r.size() > 0) rdim = *max_element(r.begin(), r.end())+1;
        return make_pair(ldim, rdim);
    }
    
    template<class Pair>
    bool compare(Pair const & p1, Pair const & p2)
    {
        return p1.first < p2.first;
    }

    struct pos_tag_lt {
        typedef std::pair<int, unsigned int> value_type;
        inline bool operator() (value_type const& lhs, value_type const& rhs)
        {
            return (lhs.first < rhs.first);
        }
    };
}

#endif
