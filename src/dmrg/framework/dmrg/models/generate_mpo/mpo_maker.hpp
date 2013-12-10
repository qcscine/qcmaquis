/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#ifndef GENERATE_MPO_MPO_MAKER_H
#define GENERATE_MPO_MPO_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/chem_compression.h"

#include "dmrg/models/lattice.h"

#include <string>
#include <sstream>

namespace generate_mpo
{
    template<class Matrix, class SymmGroup>
    class MPOMaker
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        typedef boost::tuple<size_t, size_t, op_t> block;
        typedef vector<
        pair<
        block_matrix<Matrix, SymmGroup>,
        block_matrix<Matrix, SymmGroup>
        >
        > op_pairs;
        
    public:
        MPOMaker(Lattice const& lat_,
                 const std::vector<block_matrix<Matrix, SymmGroup> > & ident_,
                 const std::vector<block_matrix<Matrix, SymmGroup> > & fill_)
        : lat(lat_)
        , length(lat.size())
        , used_dims(length)
        , identities(ident_)
        , fillings(fill_)
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
        
        void add_term(Operator_Term<Matrix, SymmGroup> const & term)
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

        MPO<Matrix, SymmGroup> create_compressed_mpo(Index<SymmGroup> const & phys, double cutoff)
        {
            if (!finalized) finalize();
            MPO<Matrix, SymmGroup> r(length);
            for (size_t p = 1; p < length - 1; ++p)
                r[p] = as_bulk(prempo[p]);
            r[0] = as_left(prempo[0]);
            r[length-1] = as_right(prempo[length-1]);

            charge_sort(get_prempo(), r);
            MPO<Matrix, SymmGroup> mpo_sorted = create_mpo();

            compressor<Matrix, SymmGroup> cpor(phys);
            MPO<Matrix, SymmGroup> mpo_out(length);
            cpor.compress(mpo_sorted, mpo_out, cutoff);

            return mpo_out;
        }

        std::vector<std::vector<block> > & get_prempo()
        {
            return prempo;
        }
        
    private:
        Lattice const& lat;
        std::size_t length;
        std::vector<block_matrix<Matrix, SymmGroup> > identities, fillings;
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
