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
	typedef block_matrix<Matrix, SymmGroup> op_t;
	typedef std::pair<op_t, bool> op_t_type;
	typedef std::pair<pos_t, op_t_type> pos_op_pair;

public:
	NTermsMPO(Lattice const& lat_,
				const std::vector<op_t> & ident_,
				const std::vector<op_t> & fill_,
				std::vector<pos_op_pair> const & op_string,
				const int & phase_)
	: lat(lat_)
	, phase(phase_)
	, identities(ident_)
	, fillings(fill_)
	{
		bool trivial_fill = true;
		for (pos_t p = 0; p < lat.size(); ++p){
			for (pos_t p_op = 0; p_op < op_string.size(); ++p_op) {
				if (p == op_string[p_op].first) {
					maquis::cout << "Inserting non-trivial operator at site " << p << std::endl;
					trivial_fill = (trivial_fill != !op_string[p_op].second.second);
					if (!trivial_fill) {
						op_t tmp;
						gemm(fillings[lat.get_prop<int>("irrep",p)], op_string[p_op].second.first, tmp);
						prempo.insert( std::make_pair(p, tmp));
					} else {
						prempo.insert( std::make_pair(p, op_string[p_op].second.first) );
					}
				}
			}
			if (prempo.count(p)==0) {
				maquis::cout << "Inserting fill/identity operator at site " << p << std::endl;
				if (trivial_fill) {
					prempo.insert( std::make_pair(p, identities[lat.get_prop<int>("irrep",p)]) );
				} else {
					prempo.insert( std::make_pair(p, fillings[lat.get_prop<int>("irrep",p)]) );
				}
			}
		}
	}

	MPO<Matrix, SymmGroup> create_mpo()
	{
		MPO<Matrix, SymmGroup> ret(prempo.size());
		MPOTensor<Matrix, SymmGroup> op(1,1);
		for (pos_t p = 0; p < prempo.size(); ++p) {
			if (phase < 0) {
				op.set(0,0, prempo[p], phase*1.0);
				phase = 1;
			} else {
				op.set(0,0, prempo[p], 1.0);
			}
			ret[p] = op;
			maquis::cout << op.at(0,0).op;
		}
		return ret;
	}

private:
	Lattice const& lat;
	int phase;
	std::vector<op_t> identities, fillings;
	std::map<pos_t, op_t> prempo;
};

