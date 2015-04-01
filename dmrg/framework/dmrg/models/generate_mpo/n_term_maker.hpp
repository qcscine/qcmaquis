/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <iostream>
#include <vector>
#include <map>

#include <alps/numeric/matrix.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

template<class Matrix, class SymmGroup>
class NTermsMPO
{
	typedef Lattice::pos_t pos_t;
	typedef typename operator_selector<Matrix, SymmGroup>::type op_t;
	typedef std::pair<pos_t, op_t> pos_op_t;

public:
	NTermsMPO(Lattice const& lat_,
				 const std::vector<op_t> & ident_,
				 const std::vector<op_t> & fill_,
				 std::vector<pos_op_t> const & ops_,
				 const int & phase)
	: lat(lat_)
	, identities(ident_)
	, fillings(fill_)
	, ops(ops_)
	{
		// Assuming all operators are fermionic!!
		bool trivial_fill=true;
		std::sort(ops.begin(),ops.end(), boost::bind(&pos_op_t::first, _1) < 
										 boost::bind(&pos_op_t::first, _2));

		for (pos_t p = 0; p < lat.size(); ++p){
			for (typename std::vector<pos_op_t>::iterator it = ops.begin(); it != ops.end(); ++it) {
				if (p == it->first) {
					if (prempo.count(p) > 0) {
						//if operator already present, multiply with it
						op_t tmp;
						if (trivial_fill) {
							op_t tmp1;
							gemm(fillings[lat.get_prop<int>("irrep",p)], it->second, tmp1);
							gemm(tmp1, prempo[p], tmp);
						} else {
							gemm(it->second, prempo[p], tmp);
						}
						prempo.erase(p);
						prempo.insert( std::make_pair(p, tmp) );
					} else {
						//else add operator to prempo
						if (trivial_fill) {
							op_t tmp;
							gemm(fillings[lat.get_prop<int>("irrep",p)], it->second, tmp);
							prempo.insert( std::make_pair(p, tmp) );
						} else {
							prempo.insert( *it );
						}
					}
					trivial_fill = !trivial_fill;
					it = ops.erase(ops.begin());
					--it;
				}
			}
			if (prempo.count(p) == 0) {
				if (trivial_fill) {
					prempo.insert( std::make_pair(p, identities[lat.get_prop<int>("irrep",p)]) );
				} else {
					prempo.insert( std::make_pair(p, fillings[lat.get_prop<int>("irrep",p)]) );
				}
			}
		}
		// Apply the phase to the first element
		prempo.begin()->second*=phase;
	}

	MPO<Matrix, SymmGroup> create_mpo()
	{
		MPO<Matrix, SymmGroup> ret(prempo.size());
		for (pos_t p = 0; p < prempo.size(); ++p) {
			MPOTensor<Matrix, SymmGroup> op(1,1);
			op.set(0,0, prempo[p], 1.0);
			ret[p] = op;
		}
		return ret;
	}
private:
	Lattice const& lat;
	std::vector<op_t> identities, fillings;
	std::vector<pos_op_t> ops;
	std::map<pos_t, op_t> prempo;
};

