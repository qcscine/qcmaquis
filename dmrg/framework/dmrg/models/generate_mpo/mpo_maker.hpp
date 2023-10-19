/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef GENERATE_MPO_MPO_MAKER_H
#define GENERATE_MPO_MPO_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/lattice/lattice.h"
#include <string>
#include <sstream>

namespace generate_mpo
{
    template<class Matrix, class SymmGroup>
    class MPOMaker
    {
        typedef typename OPTable<Matrix, SymmGroup>::op_t op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        typedef vector< pair<op_t, op_t> > op_pairs;
        
    public:
        MPOMaker(Lattice const& lat_,
                 const std::vector<op_t> & ident_,
                 const std::vector<op_t> & fill_)
        : lat(lat_)
        , length(lat.size())
        , identities(ident_)
        , fillings(fill_)
        , used_dims(length)
        , prempo(length)
        , maximum(2)
        , leftmost_right(length)
        , finalized(false)
        {
            for (size_t p = 0; p < length; ++p)
            {
                if (p+1 < length)
                    prempo[p].push_back(boost::make_tuple(std::size_t(0), std::size_t(0), identities[lat.get_prop<int>("type",p)]));
            }
        }
        
        void add_term(OperatorTerm<Matrix, SymmGroup> const & term)
        {
            // TODO: removed const&, because of sorting (non-const operation)
            std::vector<std::pair<typename Lattice::pos_t, op_t> > ops = term.operators;
            
            std::sort(ops.begin(), ops.end(), compare<std::pair<typename Lattice::pos_t, op_t> >);
            
            vector<size_t> positions;
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator
                 it = ops.begin();
                 it != ops.end(); ++it)
                positions.push_back( it->first );
            size_t minp = *min_element(positions.begin(), positions.end());
            size_t maxp = *max_element(positions.begin(), positions.end());
            
            size_t use_b = maximum++;
            
            vector<bool> done(length, false);
            for (typename vector<pair<typename Lattice::pos_t, op_t> >::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                size_t first_use_b = (it->first == minp ? 0 : use_b);
                size_t second_use_b = (it->first == maxp ? 1 : use_b);
                assert( it->first < prempo.size() );
                if (minp != maxp) { // bond term
                    prempo[it->first].push_back(boost::make_tuple(first_use_b, second_use_b, it->second));
                    used_dims[it->first].insert(use_b);
                } else // site term
                    site_terms[it->first] += it->second;
                done[it->first] = true;
            }
            
            for (size_t p = minp; p <= maxp; ++p)
                if (!done[p]) {
                    op_t const& current_filling = (term.with_sign) ? fillings[lat.get_prop<int>("type",p)] : identities[lat.get_prop<int>("type",p)];
                    prempo[p].push_back( boost::make_tuple(use_b, use_b, current_filling));
                    used_dims[p].insert(use_b);
                    done[p] = true;
                }
            
            leftmost_right = std::min(leftmost_right, maxp);
        }
        
        MPO<Matrix, SymmGroup> create_mpo()
        {
            if (!finalized) finalize(); 
            MPO<Matrix, SymmGroup> r(length);
            for (size_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[length-1] = as_right(prempo[length-1]);
            
            return r;
        }

        std::vector<std::vector<block> > & get_prempo()
        {
            return prempo;
        }
        
    private:
        Lattice const& lat;
        std::size_t length;
        std::vector<op_t> identities, fillings;
        vector<set<size_t> > used_dims;
        vector<vector<block> > prempo;
        std::map<std::size_t, op_t> site_terms;
        size_t maximum, leftmost_right;
        bool finalized;

        void finalize()
        {
            for (typename std::map<std::size_t, op_t>::const_iterator it = site_terms.begin();
                 it != site_terms.end(); ++it)
                prempo[it->first].push_back( boost::make_tuple(0, 1, it->second) );

            for (size_t p = leftmost_right + 1; p < length; ++p)
                prempo[p].push_back( boost::make_tuple(1, 1, identities[lat.get_prop<int>("type",p)]) );

            for (typename vector<vector<block> >::iterator it = prempo.begin();
                 it + 1 != prempo.end();
                 ++it)
                compress_on_bond(*it, *(it+1));

            finalized = true;
        }
        
        MPOTensor<Matrix, SymmGroup> as_bulk(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r.set(get<0>(*it), get<1>(*it), get<2>(*it));
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_left(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(1, rcd.second);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r.set(0, get<1>(*it), get<2>(*it));
            }
            return r;
        }
        
        MPOTensor<Matrix, SymmGroup> as_right(vector<block> const & ops)
        {
            pair<size_t, size_t> rcd = rcdim(ops);
            MPOTensor<Matrix, SymmGroup> r(rcd.first, 1);
            for (typename vector<block>::const_iterator it = ops.begin();
                 it != ops.end(); ++it)
            {
                r.set(get<0>(*it), 0, get<2>(*it));
            }
            return r;
        }
    };

}

#endif
