/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_BASIS_SECTOR_ITERATOR_H
#define MAQUIS_DMRG_BASIS_SECTOR_ITERATOR_H

#include "dmrg/block_matrix/indexing.h"

#include <boost/operators.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/tuple/tuple_comparison.hpp>


template <class SymmGroup>
class basis_sector_iterator_
: public boost::forward_iterator_helper<
                                          basis_sector_iterator_<SymmGroup>
                                        , std::vector<boost::tuple<typename SymmGroup::charge, std::size_t> >
                                        , std::ptrdiff_t
                                        , std::vector<boost::tuple<typename SymmGroup::charge, std::size_t> > *
                                        , std::vector<boost::tuple<typename SymmGroup::charge, std::size_t> > &
                                       >

{
    typedef typename SymmGroup::charge charge;
    typedef std::size_t size_t;
    typedef boost::tuple<charge, size_t> local_state;
    typedef typename std::vector<local_state>::const_iterator states_iterator;

    typedef const charge& (*get0_fn_t)(const boost::tuples::cons<charge, boost::tuples::cons<size_t, boost::tuples::null_type> >&);

public:
    basis_sector_iterator_()
    : valid(false)
    { }
    
    basis_sector_iterator_(size_t L_, Index<SymmGroup> const& phys, charge initc_)
    : valid(true)
    , L(L_)
    , initc(initc_)
    , it(L, 0)
    , state(L)
    {
        getter_fn = &boost::tuples::get<0, charge, boost::tuples::cons<size_t, boost::tuples::null_type> >;
    
        for (size_t i=0; i<phys.size(); ++i)
            for (size_t j=0; j<phys[i].second; ++j)
                alllocal.push_back( local_state(phys[i].first, j) );
                
        for (size_t i=0; i<L; ++i) {
            state[i] = alllocal[it[i]];
        }
        
        if (total_charge() != initc)
            advance();
    }
    
    std::vector<local_state> const& operator*() const
    {
        return state;
    }
    
    void operator++()
    {
        advance();
    }
    
    bool operator==(basis_sector_iterator_<SymmGroup> const & rhs) const
    {
        if (valid != rhs.valid)
            return false;
        if (!valid)
            return true;
        
        return (L == rhs.L) && (initc == rhs.initc) && (alllocal == rhs.alllocal) && std::equal(state.begin(), state.end(), rhs.state.begin());
    }

private:
    
    charge total_charge() const
    {
        return std::accumulate(state.begin(), state.end(), SymmGroup::IdentityCharge,
                               boost::bind(static_cast<charge(*)(charge,charge)>(&SymmGroup::fuse), _1,  boost::bind(getter_fn, _2)) );
    }
    
    void advance()
    {
        do {
            ++it[L-1];
            for (int i=L-1; (i > 0) && (it[i] == alllocal.size()); --i) {
                it[i] = 0;
                ++it[i-1];
            }
            if ( it[0] == alllocal.size() ) {
                valid = false;
                return;
            }
            
            for (size_t i=0; i<L; ++i)
                state[i] = alllocal[it[i]];
        } while(total_charge() != initc);
    }
    
    bool valid;
    size_t L;
    charge initc;
    get0_fn_t getter_fn;
    
    std::vector<local_state> alllocal;
    std::vector<size_t> it;
    std::vector<local_state> state;
};


template <class SymmGroup>
std::pair<basis_sector_iterator_<SymmGroup>, basis_sector_iterator_<SymmGroup> >
basis_sector_iterators(size_t L, Index<SymmGroup> const& phys, typename SymmGroup::charge initc=SymmGroup::IdentityCharge)
{
    return std::make_pair(basis_sector_iterator_<SymmGroup>(L, phys, initc), basis_sector_iterator_<SymmGroup>());
}


#endif
