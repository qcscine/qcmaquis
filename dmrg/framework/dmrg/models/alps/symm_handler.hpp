/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_ALPS_MODEL_SYMM_HANDLER_H
#define APP_ALPS_MODEL_SYMM_HANDLER_H

#include <alps/parameter.h>
#include <alps/model.h>

#include "dmrg/block_matrix/indexing.h"


namespace detail {
	template <class T>
	int to_integer (alps::half_integer<T> const & qn_value)
	{
		return qn_value.get_twice(); // always works with double QN, so that spin-1/2 and spin-1 are compatible
	}
}

template <class SymmGroup>
typename SymmGroup::charge state_to_charge(alps::site_state<short> const & state, alps::SiteBasisDescriptor<short> const& b,
                                           std::map<std::string, int> const& all_conserved_qn);

template <class SymmGroup>
typename SymmGroup::charge init_charge(const alps::Parameters& parms, std::map<std::string, int> const& all_conserved_qn);


template <class SymmGroup>
class symmetric_basis_descriptor {
public:
    typedef short I;
    typedef typename SymmGroup::charge charge_t;
    typedef std::size_t size_t;
    typedef std::pair<charge_t, size_t> coord_t;
    typedef std::map<std::string, int> qn_map_type;

    symmetric_basis_descriptor() : numstates_(0) { }
    
    symmetric_basis_descriptor(alps::SiteBasisDescriptor<I> const & b, qn_map_type const& all_conserved_qn)
    {
        alps::site_basis<I> states(b);

        numstates_ = states.size();
        state_index_.resize(numstates_);
        state_offset_.resize(numstates_);
        std::vector<charge_t> state_charge(numstates_);
        
        // loop over states and create phys_
        for (size_t i=0; i<states.size(); ++i) {
            state_charge[i] = state_to_charge<SymmGroup>(states[i], b, all_conserved_qn);
            
            size_t ci = phys_.position(state_charge[i]);
            if (ci < phys_.size()) phys_[ci].second += 1;
            else                   ci = phys_.insert( std::make_pair(state_charge[i], 1) );
            
            state_offset_[i] = phys_[ci].second - 1;
        }
        
        // cache mapping state -> phys
        for (size_t i=0; i<states.size(); ++i)
            state_index_[i] = phys_.position(state_charge[i]);
    }
    
    size_t size() const { return numstates_; }
    coord_t coords(size_t i) const { return coord_t(phys_[state_index_[i]].first, state_offset_[i]); }
    charge_t charge(size_t i) const { return phys_[state_index_[i]].first; }
    size_t block_size(size_t i) const { return phys_[state_index_[i]].second; }
    Index<SymmGroup> const& phys_dim() const { return phys_; }
    
private:
    size_t numstates_;
    std::vector<size_t> state_index_;
    std::vector<size_t> state_offset_;
    Index<SymmGroup> phys_;
};



#endif
