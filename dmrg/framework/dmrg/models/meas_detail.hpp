/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

#include <sstream>
#include <algorithm>

#include "dmrg/utils/utils.hpp"
#include "utils/traits.hpp"

#include <alps/hdf5.hpp>

namespace meas_detail {
    
    // forward declaration
    template<class Matrix, class SymmGroup>
	void measure_correlation_(MPS<Matrix, SymmGroup> const & mps,
							  const Lattice & lat,
							  block_matrix<Matrix, SymmGroup> const & identity,
							  block_matrix<Matrix, SymmGroup> const & fill,
							  std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
							  std::vector<std::size_t> const & order,
							  bool is_nn,
							  std::vector<typename MPS<Matrix, SymmGroup>::scalar_type>& dc,
							  std::vector<std::string>& labels,
                              bool super_meas);
    
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
        LocalMPSMeasurement(const MPS<Matrix, SymmGroup> & mps_, const Lattice & lat_)
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
                right = contraction::overlap_mpo_right_step(mps[L-i], mps[L-i], right, ident);
                right_[L-1-i] = right;
                
                left = contraction::overlap_mpo_left_step(mps[i-1], mps[i-1], left, ident);
                left_[i] = left;
            }

        }
        
        void site_term (std::pair<block_matrix<Matrix, SymmGroup>, bool> const & op,
                        std::string const & h5name,
                        std::string const & base_path) const
        {
            std::vector<typename MPSTensor<Matrix, SymmGroup>::scalar_type> vals; vals.reserve(L);
            std::vector<std::string> labels;
            MPOTensor<Matrix, SymmGroup> temp;
            temp(0,0) = op.first;

            for (int p = 0; p < L; ++p) {
                MPSTensor<Matrix, SymmGroup> vec2 =
                contraction::site_hamil2(mps[p], left_[p], right_[p], temp);
                vals.push_back( maquis::real(mps[p].scalar_overlap(vec2)) ); // MD todo: allow complex numbers
                labels.push_back( lat.get_prop<std::string>("label", p) );
            } // should return a vector of pairs or pair of vectors (todo: 30.04.12 / Matthias scalar/value types discussion)
            
            { // should be moved out to the main loop (todo: 30.04.12 / Matthias scalar/value types discussion)
                alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
                ar << alps::make_pvp(base_path + std::string("/labels"), labels);
            }
        }
        
        void bond_term (std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                        std::string const & h5name,
                        std::string const & base_path) const
        {
            assert(ops.size() == 2);
            
            std::vector<typename MPSTensor<Matrix, SymmGroup>::scalar_type> vals; vals.reserve(L-1);
            std::vector<std::string> labels;
            MPOTensor<Matrix, SymmGroup> temp;
            
            for (int p = 0; p < L-1; ++p) {
                temp(0,0) = ops[0].first;
                Boundary<Matrix, SymmGroup> tmp_b = contraction::overlap_mpo_left_step(mps[p], mps[p], left_[p], temp);
                
                temp(0,0) = ops[1].first;
                MPSTensor<Matrix, SymmGroup> vec2 =
                contraction::site_hamil2(mps[p+1], tmp_b, right_[p+1], temp);
                vals.push_back( maquis::real(mps[p+1].scalar_overlap(vec2)) ); // MD todo: allow complex numbers
                labels.push_back( lat.get_prop<std::string>("label", p, p+1) );
            } // same here (todo: 30.04.12 / Matthias scalar/value types discussion)
            
            {
                alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
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
                       std::string base_path,
                       bool super_meas = false)
	{
		std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
        std::vector<std::string> labels;
        
        if (ops.size() == 1) {
            measure_correlation_(mps, lat, identity, fill, ops, std::vector<std::size_t>(), false, vals, labels, super_meas);
        } else {
            // TODO: optimize this, by building a special MPO (MPOMaker and CorrMaker don't support it)
            for (std::size_t p = 0; p < lat.size(); ++p)
            {
            	std::vector<Lattice::pos_t> neighs = lat.forward(p);
                vals.reserve(vals.size()+neighs.size());
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
                    
                    // C - Tim, no futur, if disagree -> Alex must fix futur stuff
                    if (!super_meas){
                        vals.push_back(expval(mps, mpo));
                    }else {
                        Index<SymmGroup> phys_i = identity.left_basis();
                        typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, phys_i);
                        MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, phys_i);
                        vals.push_back(overlap(super_mpo, mps)/nn);
                    }
                    
                    labels.push_back(lat.get_prop<std::string>("label", p, *hopto));
            	}
            }
        }
        
        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(base_path + std::string("/labels"), labels);
            if ( all_true(ops.begin(), ops.end(), boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian), 
                                                              boost::bind<block_matrix<Matrix, SymmGroup> const&>(&std::pair<block_matrix<Matrix, SymmGroup>, bool>::first, _1))
                          ) ) {
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
            } else
                ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> >(1, vals));
        }
	}

    template<class Matrix, class SymmGroup>
	void measure_local_at(MPS<Matrix, SymmGroup> const & mps,
                          const Lattice & lat,
                          block_matrix<Matrix, SymmGroup> const & identity,
                          block_matrix<Matrix, SymmGroup> const & fill,
                          std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                          std::vector<std::vector<std::size_t> > const & positions,
                          std::string const & h5name,
                          std::string base_path,
                          bool super_meas = false)
	{
        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
        
        for (std::size_t p = 0; p < positions.size(); ++p)
        {
            assert( positions[p].size() == ops.size() );
            for (std::size_t i=1; i<ops.size(); ++i)
                if (positions[p][i-1] >= positions[p][i])
                    throw std::runtime_error("measure_local_at requires i1<i2<...<in.");
            
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            
            bool with_sign = false;
            for (std::size_t i=0; i<ops.size(); ++i) {
                std::size_t pos = positions[p][i];
                typedef block_matrix<Matrix, SymmGroup> op_t;
                op_t tmp;
                if (!with_sign && ops[i].second) gemm(fill, ops[i].first, tmp);
                else                             tmp = ops[i].first;
                term.operators.push_back( std::make_pair(pos, tmp) );
                
                pos++;
                with_sign = (ops[i].second) ? !with_sign : with_sign;
                if (i != ops.size()-1)
                    for (; pos<positions[p][i+1]; ++pos) {
                        if (with_sign)
                            term.operators.push_back( std::make_pair(pos, fill) );
                        else
                            term.operators.push_back( std::make_pair(pos, identity) );
                    }
            }
            term.fill_operator = identity;
            mpom.add_term(term);
            MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
            
            if (!super_meas){
                typename MPS<Matrix, SymmGroup>::scalar_type val = expval(mps, mpo);
                vals.push_back(val);
            } else {
                Index<SymmGroup> phys_i = identity.left_basis();
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, phys_i);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, phys_i);
                typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpo, mps);
                vals.push_back(val/nn);
            }
            
        }
        std::vector<std::string> labels = label_strings(lat, positions);
        
        
        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(base_path + std::string("/labels"), labels);
            if ( all_true(ops.begin(), ops.end(), boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian),
                                                              boost::bind<block_matrix<Matrix, SymmGroup> const&>(&std::pair<block_matrix<Matrix, SymmGroup>, bool>::first, _1))
                          ) ) {
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
            } else
                ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> >(1, vals));
        }
    }

    template<class Matrix, class SymmGroup>
    void measure_custom(MPS<Matrix, SymmGroup> const & mps,
                       const Lattice & lat,
                       block_matrix<Matrix, SymmGroup> const & identity,
                       block_matrix<Matrix, SymmGroup> const & fill,
                       std::vector< std::vector< std::pair<int, block_matrix<Matrix, SymmGroup> > > > const & ops,
                       std::string const & h5name,
                       std::string base_path)
    {
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat.size(), identity);

        for (int k = 0; k < ops.size(); ++k) {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators = ops[k];
            term.fill_operator = fill;
            mpom.add_term(term);
        }

        MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
        double val = expval(mps, mpo);

        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<double>(1, val));
        }
    }
 
        // C - Tim, I add one layer of wrapper for the compatibility with the futur stuff
	template<class Matrix, class SymmGroup>
	void save_helper(typename MPS<Matrix, SymmGroup>::scalar_type const& val,
    		         std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                         std::string const & h5name,
                         std::string base_path){
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            if ( all_true(ops.begin(), ops.end(), boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian), 
                                                              boost::bind<block_matrix<Matrix, SymmGroup> const&>(&std::pair<block_matrix<Matrix, SymmGroup>, bool>::first, _1))
                          ) ) {
                std::vector<double> tmp(1, maquis::real(val));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
            } else
                ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<typename MPS<Matrix, SymmGroup>::scalar_type>(1, val));
        }

   
	template<class Matrix, class SymmGroup>
	void measure_average(MPS<Matrix, SymmGroup> const & mps,
                         const Lattice & lat,
                         block_matrix<Matrix, SymmGroup> const & identity,
                         block_matrix<Matrix, SymmGroup> const & fill,
    		         std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & ops,
                         std::string const & h5name,
                         std::string base_path,
                         bool super_meas = false)
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
       
        // C - Tim, no futur, if disagree -> Alex must fix futur 
        if (!super_meas){
            typename MPS<Matrix, SymmGroup>::scalar_type  val = expval(mps, mpo);
            save_helper<Matrix, SymmGroup>(val, ops, h5name, base_path); 
        } else {
            Index<SymmGroup> phys_i = identity.left_basis();
            typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, phys_i);
            MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, phys_i);
            typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpo, mps) / nn;
            save_helper<Matrix, SymmGroup>(val, ops, h5name, base_path); 
        }

	}
    
	template<class Matrix, class SymmGroup>
	void measure_overlap(MPS<Matrix, SymmGroup> const & mps,
                         const std::string & bra_ckp,
                         std::string const & h5name, std::string const & base_path)
	{
        MPS<Matrix, SymmGroup> bra_mps;
        maquis::cout << "Measuring overlap with " << bra_ckp << "." << std::endl;
        {
            alps::hdf5::archive ar(bra_ckp, alps::hdf5::archive::READ);
            ar >> alps::make_pvp("/state", bra_mps);
        }
        
        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals(1);
        if (bra_mps[bra_mps.length()-1].col_dim().sum_of_sizes() == 1)
        {
            typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(bra_mps, mps);
            
            std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> tmp(1, val);
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
        } else {
            std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals = multi_overlap(bra_mps, mps);
            
            std::vector<std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> > tmp(1, vals);
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
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
							  std::vector<typename MPS<Matrix, SymmGroup>::scalar_type>& dc,
							  std::vector<std::string>& labels,
                              bool super_meas)
	{
        typedef boost::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

        for (size_t p = 0; p < lat.size()-(ops.size()-1); ++p) {
			std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
			std::vector<std::vector<std::size_t> > num_labels;
            maker_ptr dcorr;
            if (is_nn)
                dcorr.reset(new generate_mpo::CorrMakerNN<Matrix, SymmGroup>(mps.length(), identity, fill, ops, p) );
            else
                dcorr.reset(new generate_mpo::CorrMaker<Matrix, SymmGroup>(mps.length(), identity, fill, ops, p) );
            
            MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
            
//            maquis::cout << "site " << p << ":" << std::endl << dcorr->description() << std::endl;

            if (!super_meas)
                dct = multi_expval(mps, mpo);
            else {
                Index<SymmGroup> phys_i = identity.left_basis();
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, phys_i);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, phys_i);
                dct = multi_overlap(super_mpo, mps);
                for (int i=0; i<dct.size(); ++i)
                    dct[i] /= nn;
            }
            
            num_labels = dcorr->numeric_labels();
			std::copy(dct.begin(), dct.end(), std::back_inserter(dc));

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
							 bool is_nn=false,
                             bool super_meas=false)
	{
	    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dc;
	    std::vector<std::string> labels;
	    if (half) {
	    	measure_correlation_(mps, lat, identity, fill, ops, std::vector<std::size_t>(), is_nn, dc, labels, super_meas);
	    } else {
	    	CorrPermutator<Matrix, SymmGroup> perm(ops, is_nn);
	    	for (int i=0; i<perm.size(); ++i) {
		    	measure_correlation_(mps, lat, identity, fill, perm[i], perm.order(i), is_nn, dc, labels, super_meas);
	    	}
	    }

        {
            alps::hdf5::archive ar(h5name, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
            if ( all_true(ops.begin(), ops.end(), boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian), 
                                                              boost::bind<block_matrix<Matrix, SymmGroup> const&>(&std::pair<block_matrix<Matrix, SymmGroup>, bool>::first, _1))
                          ) ) {
                std::vector<std::vector<double> > tmp(1, maquis::real(dc));
                ar << alps::make_pvp(base_path + std::string("/mean/value"), tmp);
            } else
                ar << alps::make_pvp(base_path + std::string("/mean/value"), std::vector<std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> >(1, dc));           
            ar << alps::make_pvp(base_path + std::string("/labels"), labels);
        }
	}
    
    
} // namespace
