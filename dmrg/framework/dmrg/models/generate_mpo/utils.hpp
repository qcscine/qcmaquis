/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_UTILS_H
#define GENERATE_MPO_UTILS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

#include <boost/bind.hpp>

namespace generate_mpo
{
	template<class Matrix, class SymmGroup>
	struct Operator_Tag_Term
	{
		typedef typename OPTable<Matrix, SymmGroup>::tag_type tag_type;
        typedef typename Lattice::pos_t pos_t;
		typedef std::pair<pos_t, tag_type> op_pair_t;
        
		std::vector<op_pair_t> operators;
		tag_type fill_operator;
        typename Matrix::value_type scale;
        bool with_sign;
        
        Operator_Tag_Term() : scale(1.), with_sign(false) {}
	};
    
	template<class Matrix, class SymmGroup>
	struct Operator_Term
	{
		typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef typename Lattice::pos_t pos_t;
		typedef std::pair<pos_t, op_t> op_pair_t;
        
		std::vector<op_pair_t> operators;
		op_t fill_operator;
        bool with_sign;
        
        Operator_Term() : with_sign(false) {}
        
        void canonical_order() // TODO: check and fix for fermions
        {
            std::sort(operators.begin(), operators.end(),
                      boost::bind(&op_pair_t::first, _1) <
                      boost::bind(&op_pair_t::first, _2));
        }
        
        bool operator< (Operator_Term const & rhs) const
        {
            if (operators[0].first == rhs.operators[0].first)
                return operators.size() >= rhs.operators.size();
            return operators[0].first < rhs.operators[0].first;
        }

        bool site_match (Operator_Term const & rhs) const
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
        
        bool overlap (Operator_Term const & rhs) const
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
                return k;
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
        
        return make_pair(*max_element(l.begin(), l.end())+1,
                         *max_element(r.begin(), r.end())+1);
    }
    
    template<class Pair>
    bool compare(Pair const & p1, Pair const & p2)
    {
        return p1.first < p2.first;
    }
    
}

#endif
