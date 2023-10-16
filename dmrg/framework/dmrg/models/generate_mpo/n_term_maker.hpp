/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
							gemm(fillings[lat.get_prop<int>("type",p)], it->second, tmp1);
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
							gemm(fillings[lat.get_prop<int>("type",p)], it->second, tmp);
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
					prempo.insert( std::make_pair(p, identities[lat.get_prop<int>("type",p)]) );
				} else {
					prempo.insert( std::make_pair(p, fillings[lat.get_prop<int>("type",p)]) );
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

