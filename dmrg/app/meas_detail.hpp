/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "mp_tensors/mpo.h"
#include "mp_tensors/mps_mpo_ops.h"

#include <sstream>
#include <algorithm>

#include "utils/utils.hpp"

#include <alps/hdf5.hpp>

namespace meas_detail {
    
	inline std::vector<std::string> label_strings (const Lattice& lat, const std::vector<std::vector<std::size_t> >& labels)
	{
		std::vector<std::string> ret(labels.size());
		std::vector<std::string>::iterator ot = ret.begin();
		for (std::vector<std::vector<std::size_t> >::const_iterator it = labels.begin();
			 it != labels.end(); ++it)
		{
			std::ostringstream oss;
			for (std::vector<std::size_t>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				oss << lat.get_prop<std::string>("label", *it2);
				if (it2 + 1 != it->end())
					oss << " -- ";
			}
			*(ot++) = oss.str();
		}
		assert( ot == ret.end() );
		return ret;
	}
	inline std::vector<std::vector<std::size_t> > resort_labels (const std::vector<std::vector<std::size_t> >& labels,
																 std::vector<size_t> const & order,
																 bool is_nn=false)
	{
		std::vector<std::vector<std::size_t> > ret(labels.size());
		for (int i=0; i<labels.size(); ++i) {
#ifndef NDEBUG
			if (is_nn)
				assert(2*order.size() == labels[i].size());
			else
				assert(order.size() == labels[i].size());
#endif
			ret[i].resize(order.size());


			for (int j=0; j<order.size(); ++j) {
				if (is_nn) {
					ret[i][order[j]] = labels[i][2*j];
					ret[i][order[j]+1] = labels[i][2*j+1];
				} else {
					ret[i][order[j]] = labels[i][j];
				}
			}
		}
		return ret;
	}
    
	template<class Matrix, class SymmGroup>
	void measure_local(MPS<Matrix, SymmGroup> & mps,
                       const Lattice & lat,
                       block_matrix<Matrix, SymmGroup> const & identity,
                       block_matrix<Matrix, SymmGroup> const & fill,
						  std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                       alps::hdf5::oarchive & ar,
                       std::string base_path)
	{
		std::vector<double> vals;
        std::vector<std::string> labels;
		for (std::size_t p = 0; p < lat.size(); ++p)
        {
            if (ops.size() == 1) {
				generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
				generate_mpo::Operator_Term<Matrix, SymmGroup> term;
				term.operators.push_back( std::make_pair(p, ops[0].first) );
				mpom.add_term(term);
				MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();

				double val = expval(mps, mpo, 0);
				vals.push_back(val);
				labels.push_back(lat.get_prop<std::string>("label", p));
            } else {
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
            	for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
            		 hopto != neighs.end();
            		 ++hopto)
            	{
					generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
					generate_mpo::Operator_Term<Matrix, SymmGroup> term;
					term.operators.push_back( std::make_pair(p, ops[0].first) );
					term.operators.push_back( std::make_pair(*hopto, ops[1].first) );
					term.fill_operator = (ops[0].second) ? fill : identity;
					mpom.add_term(term);
					MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();

					double val = expval(mps, mpo, 0);
					vals.push_back(val);
					labels.push_back(lat.get_prop<std::string>("label", p, *hopto));
            	}
            }
        }
        
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), vals);
	    ar << alps::make_pvp(base_path + std::string("/labels"), labels);
	}
    
	template<class Matrix, class SymmGroup>
	void measure_average(MPS<Matrix, SymmGroup> & mps,
                         const Lattice & lat,
                         block_matrix<Matrix, SymmGroup> const & identity,
                         block_matrix<Matrix, SymmGroup> const & fill,
						 std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                         alps::hdf5::oarchive & ar,
                         std::string base_path)
	{
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
		for (std::size_t p = 0; p < lat.size(); ++p)
        {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators.push_back( std::make_pair(p, ops[0].first) );
            if (ops.size() == 1) {
                mpom.add_term(term);
            } else {
				term.fill_operator = (ops[0].second) ? fill : identity;
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
            	for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
            		 hopto != neighs.end();
            		 ++hopto)
            	{
                    generate_mpo::Operator_Term<Matrix, SymmGroup> term2(term);
                    term2.operators.push_back( std::make_pair(*hopto, ops[1].first) );
            	}

            }
        }
        MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
        double val = expval(mps, mpo, 0);
        
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), val);
	}

	template <class Matrix, class SymmGroup>
	class CorrPermutator
	{
	public:
		typedef std::size_t size_t;
		typedef std::pair<block_matrix<Matrix, SymmGroup>, bool> inner_t;
		typedef std::vector<inner_t> value_t;
		CorrPermutator (value_t const & ops, bool is_nn)
		{
			std::vector<size_t> ord;
			if (is_nn)
				for (size_t i=0; i<ops.size(); i+=2) ord.push_back(i);
			else
				for (size_t i=0; i<ops.size(); ++i) ord.push_back(i);

			do {
				std::vector<size_t> check_pre;
				for (size_t i=0; i<ops.size(); ++i) {
					if (ops[i].second) {
						if (!is_nn)
							check_pre.push_back(ord[i]);
						else
							if (i % 2 == 0)
								check_pre.push_back(ord[i]);
							else
								check_pre.push_back(ord[i]+1);
					}
				}

				cmp_with_prefactor cmp;
				cmp.prefactor = 1.;
				std::sort(check_pre.begin(), check_pre.end(), cmp);

				value_t tmp;
				for (size_t i=0; i<ord.size(); ++i) {
					tmp.push_back( ops[ord[i]] );
					if (is_nn) {
						tmp.push_back( ops[ord[i]+1] );
					}
				}
				tmp[0].first *= cmp.prefactor;
				perm.push_back(tmp);
				orders.push_back(ord);
			} while (std::next_permutation(ord.begin(), ord.end()));
		}

		size_t size() const {return perm.size();}
		value_t operator[](size_t i) const {return perm[i];}
		std::vector<size_t> order(size_t i) const {return orders[i];}

	private:
		std::vector<value_t> perm;
		std::vector<std::vector<size_t> > orders;
	};

	template<class Matrix, class SymmGroup>
	void measure_correlation_(MPS<Matrix, SymmGroup> & mps,
							  const Lattice & lat,
							  block_matrix<Matrix, SymmGroup> const & identity,
							  block_matrix<Matrix, SymmGroup> const & fill,
							  std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
							  std::vector<std::size_t> const & order,
							  bool is_nn,
							  std::vector<double>& dc,
							  std::vector<std::string>& labels)
	{
		for (size_t p = 0; p < lat.size()-(ops.size()-1); ++p) {
			std::vector<double> dct;
			std::vector<std::vector<std::size_t> > num_labels;
			if (is_nn) {
				generate_mpo::CorrMakerNN<Matrix, SymmGroup> dcorr(mps.length(), identity, fill,
																   ops, p);
				MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
				dct = multi_expval(mps, mpo);
				num_labels = dcorr.numeric_labels();
			} else {
				generate_mpo::CorrMaker<Matrix, SymmGroup> dcorr(mps.length(), identity, fill,
																 ops, p);
				MPO<Matrix, SymmGroup> mpo = dcorr.create_mpo();
				dct = multi_expval(mps, mpo);
				num_labels = dcorr.numeric_labels();
			}
			std::copy(dct.begin(), dct.end(), std::back_inserter(dc));

			//std::cout << dcorr.description() << std::endl;


			std::vector<std::string> lbt;
			if (order.size() > 0)
				lbt = label_strings(lat, resort_labels(num_labels, order, is_nn));
			else
				lbt = label_strings(lat, num_labels);
			std::copy(lbt.begin(), lbt.end(), std::back_inserter(labels));
		}
	}

	template<class Matrix, class SymmGroup>
	void measure_correlation(MPS<Matrix, SymmGroup> & mps,
							 const Lattice & lat,
							 block_matrix<Matrix, SymmGroup> const & identity,
							 block_matrix<Matrix, SymmGroup> const & fill,
							 std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
							 alps::hdf5::oarchive & ar,
							 std::string base_path,
							 bool half=false,
							 bool is_nn=false)
	{
	    std::vector<double> dc;
	    std::vector<std::string> labels;
	    if (half) {
	    	measure_correlation_(mps, lat, identity, fill, ops, std::vector<std::size_t>(), is_nn, dc, labels);
	    } else {
	    	CorrPermutator<Matrix, SymmGroup> perm(ops, is_nn);
	    	for (int i=0; i<perm.size(); ++i) {
		    	measure_correlation_(mps, lat, identity, fill, perm[i], perm.order(i), is_nn, dc, labels);
	    	}
	    }
	    ar << alps::make_pvp(base_path + std::string("/mean/value"), dc);
	    ar << alps::make_pvp(base_path + std::string("/labels"), labels);
	}
    
    
} // namespace
