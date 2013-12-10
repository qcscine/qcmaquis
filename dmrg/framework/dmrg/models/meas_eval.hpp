/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"
#include "dmrg/models/meas_prepare.hpp"
#include "dmrg/models/generate_mpo.hpp"

#include <sstream>
#include <algorithm>
#include <boost/iterator/counting_iterator.hpp>

#include "dmrg/utils/utils.hpp"
#include "utils/traits.hpp"

namespace meas_eval {
    
	inline std::vector<std::string> label_strings (const Lattice& lat, const std::vector<std::vector<std::size_t> >& labels)
	{
        std::vector<std::string> ret;
        ret.reserve(labels.size());
		for (std::vector<std::vector<std::size_t> >::const_iterator it = labels.begin();
			 it != labels.end(); ++it)
		{
			std::ostringstream oss;
			for (std::vector<std::size_t>::const_iterator it2 = it->begin(); it2 != it->end(); ++it2) {
				oss << lat.get_prop<std::string>("label", *it2);
				if (it2 + 1 != it->end())
					oss << " -- ";
			}
            ret.push_back(oss.str());
		}
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
			ret[i].resize(labels[i].size());


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
    bool is_hermitian_meas(std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops)
    {
        bool is_herm = true;
        for (int i=0; i<ops.size() && is_herm; ++i)
            is_herm = all_true(ops[i].first.begin(), ops[i].first.end(),
                               boost::bind(static_cast<bool (*)(block_matrix<Matrix, SymmGroup> const&)>(&is_hermitian), _1));
        return is_herm;
    }
    
    template<class T>
	void save_helper(T const& val, bool cast_to_real,
                     std::string const & h5name, std::string const& base_path, std::string const& obs_name)
    {
        storage::archive ar(h5name, "w");
        
        if (cast_to_real)
        {
            std::vector<double> tmp(1, maquis::real(val));
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
        } else {
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << std::vector<T>(1, val);
        }
    }

    template<class T>
	void save_helper(std::vector<T> const& val, std::vector<std::string> const& labels, bool cast_to_real,
                     std::string const & h5name, std::string const& base_path, std::string const& obs_name)
    {
        storage::archive ar(h5name, "w");
        
        ar[base_path + storage::encode(obs_name) + std::string("/labels")] << labels;
        if (cast_to_real)
        {
            std::vector<std::vector<double> > tmp(1, maquis::real(val));
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
        } else {
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << std::vector<std::vector<T> >(1, val);
        }
    }

    
    template <class Matrix, class SymmGroup>
    class LocalMPSMeasurement
    {
        typedef typename SymmGroup::subcharge subcharge;
    public:
        LocalMPSMeasurement(const MPS<Matrix, SymmGroup> & mps_, const Lattice & lat_)
        : mps(mps_)
        , lat(lat_)
        , L(mps.size())
        , left_(L)
        , right_(L)
        , initialized(false)
        { }
        
        void init() const
        {
            // init right_ & left_
            Boundary<Matrix, SymmGroup> right = mps.right_boundary(), left = mps.left_boundary();
            right_[L-1] = right;
            left_[0] = left;
            for (int i = 1; i < L; ++i) {
                {
                    MPOTensor<Matrix, SymmGroup> ident;
                    ident.set(0, 0, identity_matrix<Matrix>(mps[L-i].site_dim()));
                    right = contraction::overlap_mpo_right_step(mps[L-i], mps[L-i], right, ident);
                    right_[L-1-i] = right;
                }
                {
                    MPOTensor<Matrix, SymmGroup> ident;
                    ident.set(0, 0, identity_matrix<Matrix>(mps[i-1].site_dim()));
                    left = contraction::overlap_mpo_left_step(mps[i-1], mps[i-1], left, ident);
                    left_[i] = left;
                }
            }
        }
        
        void site_term (std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> const & op,
                        std::string const & h5name,
                        std::string const & base_path,
                        std::string const & obs_name) const
        {
            if (!initialized)
                init();
            
            std::vector<typename MPSTensor<Matrix, SymmGroup>::scalar_type> vals; vals.reserve(L);
            std::vector<std::string> labels;

            for (typename Lattice::pos_t p = 0; p < L; ++p) {
                subcharge type = lat.get_prop<subcharge>("type", p);
                if (op.first[type].n_blocks() > 0) {
                    MPOTensor<Matrix, SymmGroup> temp;
                    temp.set(0, 0, op.first[type]);
                    
                    MPSTensor<Matrix, SymmGroup> vec2 =
                    contraction::site_hamil2(mps[p], left_[p], right_[p], temp);
                    vals.push_back( maquis::real(mps[p].scalar_overlap(vec2)) ); // MD todo: allow complex numbers
                    labels.push_back( lat.get_prop<std::string>("label", p) );
                }
            } // should return a vector of pairs or pair of vectors (todo: 30.04.12 / Matthias scalar/value types discussion)
            
            { // should be moved out to the main loop (todo: 30.04.12 / Matthias scalar/value types discussion)
                storage::archive ar(h5name, "w");
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
                ar[base_path + storage::encode(obs_name) + std::string("/labels")] << labels;
            }
        }
        
        void bond_term (std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops,
                        std::string const & h5name,
                        std::string const & base_path,
                        std::string const & obs_name) const
        {
            assert(ops.size() == 2);
            
            if (!initialized)
                init();
            
            std::vector<typename MPSTensor<Matrix, SymmGroup>::scalar_type> vals; vals.reserve(L-1);
            std::vector<std::string> labels;
            MPOTensor<Matrix, SymmGroup> temp;
            
            for (int p = 0; p < L-1; ++p) {
                int type1 = lat.get_prop<int>("type",p);
                int type2 = lat.get_prop<int>("type",p);
                if (ops[0].first[type1].n_blocks() > 0 && ops[1].first[type2].n_blocks() > 0) {
                    temp.set(0, 0, ops[0].first[type1]);
                    Boundary<Matrix, SymmGroup> tmp_b = contraction::overlap_mpo_left_step(mps[p], mps[p], left_[p], temp);
                    
                    temp.set(0, 0, ops[1].first[type2]);
                    MPSTensor<Matrix, SymmGroup> vec2 =
                    contraction::site_hamil2(mps[p+1], tmp_b, right_[p+1], temp);
                    vals.push_back( maquis::real(mps[p+1].scalar_overlap(vec2)) ); // MD todo: allow complex numbers
                    labels.push_back( lat.get_prop<std::string>("label", p, p+1) );
                }
            } // same here (todo: 30.04.12 / Matthias scalar/value types discussion)
            
            {
                storage::archive ar(h5name, "w");
                std::vector<std::vector<double> > tmp(1, maquis::real(vals));
                ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
                ar[base_path + storage::encode(obs_name) + std::string("/labels")] << labels;
            }
        }
        
    private:
        const MPS<Matrix, SymmGroup> & mps;
        const Lattice & lat;
        int L;
        mutable std::vector<Boundary<Matrix, SymmGroup> > left_, right_;
        bool initialized;

    };
    
    template<class Matrix, class SymmGroup>
	void measure_local(MPS<Matrix, SymmGroup> const & mps,
                       const Lattice & lat,
                       std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                       std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                       Measurement_Term<Matrix, SymmGroup> const& term,
                       std::string const & h5name,
                       std::string base_path,
                       bool super_meas = false)
    {
        std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops = term.operators;

        std::vector<std::string> labels;
        std::vector<MPO<Matrix, SymmGroup> > mpos;
        
        boost::tie(mpos, labels) = meas_prepare::local(lat, identities, fillings, ops);
        
        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals; vals.reserve(mpos.size());
        typename MPS<Matrix, SymmGroup>::scalar_type nn;

        if (super_meas)
            nn = dm_trace(mps, term.phys_psi);
            
        for (std::size_t i=0; i<mpos.size(); ++i) {
            if (!super_meas) {
                vals.push_back( expval(mps, mpos[i]) );
            } else {
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpos[i], term.phys_psi);
                typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpo, mps);
                vals.push_back(val/nn);
            }
        }
        
        save_helper(vals, labels, is_hermitian_meas(term.operators), h5name, base_path, term.name);
	}

    template<class Matrix, class SymmGroup>
	void measure_local_at(MPS<Matrix, SymmGroup> const & mps,
                          const Lattice & lat,
                          std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                          std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                          Measurement_Term<Matrix, SymmGroup> const& term,
                          std::string const & h5name,
                          std::string base_path,
                          bool super_meas = false)
	{
        std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> > , bool> > const & ops = term.operators;
        std::vector<std::vector<std::size_t> > const & positions = term.positions;

        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
        
        for (std::size_t p = 0; p < positions.size(); ++p)
        {
            assert( positions[p].size() == ops.size() );
            for (std::size_t i=1; i<ops.size(); ++i)
                if (positions[p][i-1] >= positions[p][i])
                    throw std::runtime_error("measure_local_at requires i1<i2<...<in.");
            
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
            generate_mpo::Operator_Term<Matrix, SymmGroup> hterm;
            
            bool with_sign = false;
            for (std::size_t i=0; i<ops.size(); ++i) {
                std::size_t pos = positions[p][i];
                
                block_matrix<Matrix, SymmGroup> const& fill  = fillings[lat.get_prop<int>("type", pos)];
                block_matrix<Matrix, SymmGroup> const& op    = ops[i].first[lat.get_prop<int>("type", pos)];
                
                typedef block_matrix<Matrix, SymmGroup> op_t;
                op_t tmp;
                if (!with_sign && ops[i].second) gemm(fill, op, tmp);
                else                             tmp = op;
                hterm.operators.push_back( std::make_pair(pos, tmp) );
                
                pos++;
                with_sign = (ops[i].second) ? !with_sign : with_sign;
                if (i != ops.size()-1)
                    for (; pos<positions[p][i+1]; ++pos) {
                        block_matrix<Matrix, SymmGroup> const& fill  = fillings[lat.get_prop<int>("type", pos)];
                        block_matrix<Matrix, SymmGroup> const& ident = identities[lat.get_prop<int>("type", pos)];
                        hterm.operators.push_back( std::make_pair(pos, (with_sign) ? fill : ident) );
                    }
            }
            hterm.with_sign = false; // filling already taken care above
            mpom.add_term(hterm);
            MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
            
            if (!super_meas){
                typename MPS<Matrix, SymmGroup>::scalar_type val = expval(mps, mpo);
                vals.push_back(val);
            } else {
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, term.phys_psi);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, term.phys_psi);
                typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpo, mps);
                vals.push_back(val/nn);
            }
            
        }
        std::vector<std::string> labels = label_strings(lat, positions);
        save_helper(vals, labels, is_hermitian_meas(ops), h5name, base_path, term.name);
    }

    template<class Matrix, class SymmGroup>
    void measure_custom(MPS<Matrix, SymmGroup> const & mps,
                        const Lattice & lat,
                        std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                        std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                        std::vector< std::vector< std::pair<int, block_matrix<Matrix, SymmGroup> > > > const & ops,
                        std::string const & h5name,
                        std::string base_path,
                        std::string const & obs_name)
    {
        generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);

        for (int k = 0; k < ops.size(); ++k) {
            generate_mpo::Operator_Term<Matrix, SymmGroup> term;
            term.operators = ops[k];
            term.with_sign = true;
            mpom.add_term(term);
        }

        MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
        double val = expval(mps, mpo);

        {
            storage::archive ar(h5name, "w");
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << std::vector<double>(1, val);
        }
    }
   
	template<class Matrix, class SymmGroup>
	void measure_average(MPS<Matrix, SymmGroup> const & mps,
                         const Lattice & lat,
                         std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                         std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                         Measurement_Term<Matrix, SymmGroup> const& term,
                         std::string const & h5name,
                         std::string base_path,
                         bool super_meas = false)
	{
        MPO<Matrix, SymmGroup> mpo = meas_prepare::average(lat, identities, fillings, term.operators);
       
        // C - Tim, no futur, if disagree -> Alex must fix futur 
        if (!super_meas){
            typename MPS<Matrix, SymmGroup>::scalar_type  val = expval(mps, mpo);
            save_helper(val, is_hermitian_meas(term.operators), h5name, base_path, term.name);
        } else {
            typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, term.phys_psi);
            MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, term.phys_psi);
            typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpo, mps) / nn;
            save_helper(val, is_hermitian_meas(term.operators), h5name, base_path, term.name);
        }

	}
    
    template<class Matrix, class SymmGroup>
	void dm_overlap(MPS<Matrix, SymmGroup> const & mps, MPS<Matrix, SymmGroup> const& mps_ident,
                    std::vector<MPS<Matrix, SymmGroup> > const& super_mpos,
                    std::vector<std::string> const& labels,
                    bool is_multi_overlap,
                    std::string const & h5name, std::string const & base_path, std::string const & obs_name)
	{
        typename MPS<Matrix, SymmGroup>::scalar_type nn = overlap(mps, mps_ident);
        
        if (labels.size() == 0) {
            typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(super_mpos[0], mps) / nn;
            std::vector<double> tmp(1, maquis::real(val));
            
            storage::archive ar(h5name, "w");
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
        } else {
            std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals;
            for (std::size_t i=0; i<super_mpos.size(); ++i) {
                if (is_multi_overlap) {
                    std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> tmp = multi_overlap(super_mpos[i], mps);
                    for(int i = 0; i < tmp.size(); ++i) tmp[i] /= nn;
                    std::copy(tmp.begin(), tmp.end(), std::back_inserter(vals));
                } else {
                    vals.push_back( overlap(super_mpos[i], mps) / nn );
                }
                
            }
            std::vector<std::vector<double> > tmp(1, maquis::real(vals));
            storage::archive ar(h5name, "w");
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
            ar[base_path + storage::encode(obs_name) + std::string("/labels")] << labels;
        }
	}

	template<class Matrix, class SymmGroup>
	void measure_overlap(MPS<Matrix, SymmGroup> const & mps,
                         const std::string & bra_ckp,
                         std::string const & h5name, std::string const & base_path, std::string const & obs_name)
	{
        maquis::cout << "Measuring overlap with " << bra_ckp << "." << std::endl;
        MPS<Matrix, SymmGroup> bra_mps;
        load(bra_ckp, bra_mps);
        
        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals(1);
        if (bra_mps[bra_mps.length()-1].col_dim().sum_of_sizes() == 1)
        {
            typename MPS<Matrix, SymmGroup>::scalar_type val = overlap(bra_mps, mps);
            
            std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> tmp(1, val);
            storage::archive ar(h5name, "w");
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
        } else {
            std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> vals = multi_overlap(bra_mps, mps);
            
            std::vector<std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> > tmp(1, vals);
            storage::archive ar(h5name, "w");
            ar[base_path + storage::encode(obs_name) + std::string("/mean/value")] << tmp;
        }
	}
    
	template <class Matrix, class SymmGroup>
	class CorrPermutator
	{
	public:
		typedef std::size_t size_t;
		typedef std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> inner_t;
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
				
                for (int type=0; type<tmp[0].first.size(); ++type)
                    tmp[0].first[type] *= cmp.prefactor;
                
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
                              std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                              std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                              std::vector<std::size_t> const& positions_first,
							  std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops,
							  std::vector<std::size_t> const & order,
							  bool is_nn,
							  std::vector<typename MPS<Matrix, SymmGroup>::scalar_type>& dc,
							  std::vector<std::string>& labels,
                              Index<SymmGroup> const& phys_psi,
                              bool super_meas)
	{
        typedef boost::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

//        for (size_t p = 0; p < lat.size()-(ops.size()-1); ++p) {
        for (std::vector<std::size_t>::const_iterator it = positions_first.begin();
             it != positions_first.end(); ++it) {
            if (*it >= lat.size()-(ops.size()-1))
                throw std::runtime_error("cannot measure correlation with first operator at p="+boost::lexical_cast<std::string>(*it)+".");
            maquis::cout << "  site " << *it << std::endl;
			std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct;
			std::vector<std::vector<std::size_t> > num_labels;
            maker_ptr dcorr;
            if (is_nn)
                dcorr.reset(new generate_mpo::CorrMakerNN<Matrix, SymmGroup>(lat, identities, fillings, ops, *it) );
            else
                dcorr.reset(new generate_mpo::CorrMaker<Matrix, SymmGroup>(lat, identities, fillings, ops, *it) );
            
            MPO<Matrix, SymmGroup> mpo = dcorr->create_mpo();
            
//            maquis::cout << "site " << p << ":" << std::endl << dcorr->description() << std::endl;

            if (!super_meas)
                dct = multi_expval(mps, mpo);
            else {
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, phys_psi);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, phys_psi);
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
                             std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
                             std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
                             Measurement_Term<Matrix, SymmGroup> const& term,
							 std::string const & h5name,
							 std::string base_path,
							 bool half=false,
							 bool is_nn=false,
                             bool super_meas=false)
	{
        std::vector<std::vector<std::size_t> > const & positions = term.positions;
        std::vector<std::pair<std::vector<block_matrix<Matrix, SymmGroup> >, bool> > const & ops = term.operators;

        std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dc;
	    std::vector<std::string> labels;
        
        assert(positions.size() < 2);
        std::vector<std::size_t> positions_first;
        if (positions.size() == 1) {
            positions_first = positions[0];
        } else {
            std::copy(boost::counting_iterator<int>(0), boost::counting_iterator<int>(lat.size()-(ops.size()-1)),
                      back_inserter(positions_first));
        }
        
	    if (half) {
	    	measure_correlation_(mps, lat, identities, fillings, positions_first, ops, std::vector<std::size_t>(), is_nn, dc, labels, term.phys_psi, super_meas);
	    } else {
	    	CorrPermutator<Matrix, SymmGroup> perm(ops, is_nn);
	    	for (int i=0; i<perm.size(); ++i) {
		    	measure_correlation_(mps, lat, identities, fillings, positions_first, perm[i], perm.order(i), is_nn, dc, labels, term.phys_psi, super_meas);
	    	}
	    }
        
        save_helper(dc, labels, is_hermitian_meas(ops), h5name, base_path, term.name);
	}
    
    
} // namespace
