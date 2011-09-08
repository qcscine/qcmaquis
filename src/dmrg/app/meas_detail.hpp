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
#include <alps/numeric/real.hpp>

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
    
    template <class Matrix, class SymmGroup>
    class LocalMPSMeasurement
    {
    public:
        LocalMPSMeasurement (const MPS<Matrix, SymmGroup> & mps_, const Lattice & lat_)
        : mps(mps_)
        , lat(lat_)
        , L(mps.size())
        , phys_i(mps[0].site_dim())
        , left_(L)
        , right_(L)
        {
            ident(0,0) = identity_matrix<Matrix>(phys_i);
            
            // init right_ & left_
            Boundary<Matrix, SymmGroup> right = mps.right_boundary(), left = mps.left_boundary();
            right_[L-1] = right;
            left_[0] = left;
            for (int i = 1; i < L; ++i) {
                MPSTensor<Matrix, SymmGroup> bkp;

                bkp = mps[L-i];
                right = contraction::overlap_mpo_right_step(mps[L-i], bkp, right, ident);
                right_[L-1-i] = right;
                
                bkp = mps[i-1];
                left = contraction::overlap_mpo_left_step(mps[i-1], bkp, left, ident);
                left_[i] = left;
            }

        }
        
        void site_term (std::pair<block_matrix<Matrix, SymmGroup>, bool> const & op,
                        std::string const & h5name,
                        std::string const & base_path) const
        {
            std::vector<double> vals;
            std::vector<std::string> labels;
            MPOTensor<Matrix, SymmGroup> temp;

            for (int p = 0; p < L; ++p) {
                temp(0,0) = op.first;
                MPSTensor<Matrix, SymmGroup> vec2 =
                contraction::site_hamil2(mps[p], left_[p], right_[p], temp);
                vals.push_back( alps::numeric::real(mps[p].scalar_overlap(vec2)) );
                labels.push_back( lat.get_prop<std::string>("label", p) );
            }
            
            {
                alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE);
                ar << alps::make_pvp(base_path + std::string("/mean/value"), vals);
                ar << alps::make_pvp(base_path + std::string("/labels"), labels);
            }
        }
        
        void bond_term (std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                        std::string const & h5name,
                        std::string const & base_path) const
        {
            assert(ops.size() == 2);
            
            std::vector<double> vals;
            std::vector<std::string> labels;
            MPOTensor<Matrix, SymmGroup> temp;
            Boundary<Matrix, SymmGroup> tmp_b;
            
            for (int p = 0; p < L-1; ++p) {
                temp(0,0) = ops[0].first;
                MPSTensor<Matrix, SymmGroup> bkp = mps[p];
                
                tmp_b = contraction::overlap_mpo_left_step(mps[p], bkp, left_[p], temp);
                
                temp(0,0) = ops[1].first;
                MPSTensor<Matrix, SymmGroup> vec2 =
                contraction::site_hamil2(mps[p+1], tmp_b, right_[p+1], temp);
                vals.push_back( alps::numeric::real(mps[p+1].scalar_overlap(vec2)) );
                labels.push_back( lat.get_prop<std::string>("label", p, p+1) );
            }
            
            {
                alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE);
                ar << alps::make_pvp(base_path + std::string("/mean/value"), vals);
                ar << alps::make_pvp(base_path + std::string("/labels"), labels);
            }
        }
        
    private:
        const MPS<Matrix, SymmGroup> & mps;
        const Lattice & lat;
        int L;
        Index<SymmGroup> phys_i;
        MPOTensor<Matrix, SymmGroup> ident;
        std::vector<Boundary<Matrix, SymmGroup> > left_, right_;

    };
    
    template<class Matrix, class SymmGroup>
	void measure_local(MPS<Matrix, SymmGroup> const & mps,
                       const Lattice & lat,
                       block_matrix<Matrix, SymmGroup> const & identity,
                       block_matrix<Matrix, SymmGroup> const & fill,
                       std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                       std::string const & h5name,
                       std::string base_path)
	{
		std::vector<double> vals;
        std::vector<std::string> labels;
        
        if (ops.size() == 1) {
            measure_correlation_(mps, lat, identity, fill, ops, std::vector<std::size_t>(), false, vals, labels);
        } else {
            // TODO: optimize this, by building a special MPO (MPOMaker and CorrMaker don't support it)
            for (std::size_t p = 0; p < lat.size(); ++p)
            {
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
                    
					double val = expval(mps, mpo);
					vals.push_back(val);
					labels.push_back(lat.get_prop<std::string>("label", p, *hopto));
            	}
            }
        }
                
        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), vals);
            ar << alps::make_pvp(base_path + std::string("/labels"), labels);
        }
	}

	template<class Matrix, class SymmGroup>
	void measure_average(MPS<Matrix, SymmGroup> const & mps,
                         const Lattice & lat,
                         block_matrix<Matrix, SymmGroup> const & identity,
                         block_matrix<Matrix, SymmGroup> const & fill,
						 std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                         std::string const & h5name,
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
        double val = expval(mps, mpo);
        
        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<double>(1, val));
        }
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
	void measure_correlation_(MPS<Matrix, SymmGroup> const & mps,
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
	void measure_correlation(MPS<Matrix, SymmGroup> const & mps,
							 const Lattice & lat,
							 block_matrix<Matrix, SymmGroup> const & identity,
							 block_matrix<Matrix, SymmGroup> const & fill,
							 std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
							 std::string const & h5name,
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
        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), dc);
            ar << alps::make_pvp(base_path + std::string("/labels"), labels);
        }
	}
    
    
} // namespace
