/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
 *               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MEASUREMENTS_REL_NRANKRDM_H
#define MEASUREMENTS_REL_NRANKRDM_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/measurements/correlations.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"
#include "dmrg/models/generate_mpo/n_term_maker.hpp"

#include <algorithm>
#include <functional>
#include <boost/iterator/counting_iterator.hpp>
#include <boost/filesystem.hpp>

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class Rel_NRankRDM : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef typename operator_selector<Matrix, SymmGroup>::type op_t;
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        typedef Lattice::pos_t pos_t;
		typedef std::pair<op_t, bool> op_t_type;
        typedef std::vector<pos_t> positions_type;
		typedef std::pair<pos_t, op_t> pos_op_t;
    
    public:
        Rel_NRankRDM(std::string const& name_, const Lattice & lat,
                 op_vec const & identities_, op_vec const & fillings_,
                 std::vector<bond_element> const& ops_,
                 bool half_only_, bool nearest_neighbors_only,
                 positions_type const& positions_ = positions_type(),
                 std::string const& ckp_ = std::string(""))
        : base(name_)
        , lattice(lat)
        , positions_first(positions_)
        , identities(identities_)
        , fillings(fillings_)
        , ops(ops_)
        , half_only(half_only_)
        , is_nn(nearest_neighbors_only)
        , bra_ckp(ckp_)
        {
            pos_t extent = ops.size() > 2 ? lattice.size() : lattice.size()-1;
            if (positions_first.size() == 0)
                std::copy(boost::counting_iterator<pos_t>(0), boost::counting_iterator<pos_t>(extent),
                          back_inserter(positions_first));
            
            this->cast_to_real = is_hermitian_meas(ops[0]);
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& ket_mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();

            MPS<Matrix, SymmGroup> bra_mps;
            if (bra_ckp != "") {
                if(boost::filesystem::exists(bra_ckp))
                    load(bra_ckp, bra_mps);
                else
                    throw std::runtime_error("The bra checkpoint file " + bra_ckp + " was not found\n");
            }

            if (ops[0].size() == 2)
                measure_correlation(bra_mps, ket_mps);
            else if (ops[0].size() == 4)
                //measure_2rdm(bra_mps, ket_mps, ops);
                measure_2rdm(bra_mps, ket_mps);
            else
                throw std::runtime_error("rel_nrankrdm correlation measurements at the moment supported with 2 and 4 operators");
        }
        
    protected:
        typedef boost::shared_ptr<generate_mpo::CorrMakerBase<Matrix, SymmGroup> > maker_ptr;

        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new Rel_NRankRDM(*this);
        }
        
        void measure_correlation(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                                 MPS<Matrix, SymmGroup> const & ket_mps
                                )
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            //int ecounter = 0;

            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for schedule(dynamic)
            #endif
		for (pos_t i = 0; i < lattice.size(); ++i){
			for (pos_t j = i; j < lattice.size(); ++j){

				int phase = 1;
				std::vector<pos_op_t> op_string;
				op_string.push_back( std::make_pair(i, ops[0][0].first[lattice.get_prop<int>("type",i)]));
				op_string.push_back( std::make_pair(j, ops[0][1].first[lattice.get_prop<int>("type",j)]));

				////// CALL MPO MAKER /////////
				NTermsMPO<Matrix, SymmGroup> rdm_elem(lattice, identities, fillings, op_string, phase);
				MPO<Matrix, SymmGroup> mpo = rdm_elem.create_mpo();
				//////////////////////////////////

				//std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(bra_mps, ket_mps, mpo);
				typename MPS<Matrix, SymmGroup>::scalar_type dct = expval(bra_mps, ket_mps, mpo);
						
                                         //ecounter++;
                                          
				if(dct != 0.0) {
					//maquis::cout << std::fixed << std::setprecision(10) << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << "\t" << dct << std::endl;
							
							
					std::vector<pos_t> label; label.push_back(i); label.push_back(j);
					std::vector<std::vector<pos_t> > num_labels;
					num_labels.push_back(label);
					std::vector<std::string> lbt = label_strings(lattice, num_labels);
								
					this->vector_results.push_back(dct);
					this->labels.push_back(lbt[0]);
				} 
			} // j loop
		} // i loop
        }

        void measure_2rdm(MPS<Matrix, SymmGroup> const & dummy_bra_mps,
                          MPS<Matrix, SymmGroup> const & ket_mps
                          //std::vector<bond_element> const & ops,
                         )
        {
            // Test if a separate bra state has been specified
            bool bra_neq_ket = (dummy_bra_mps.length() > 0);
            MPS<Matrix, SymmGroup> const & bra_mps = (bra_neq_ket) ? dummy_bra_mps : ket_mps;

            //int ecounter = 0;

            // TODO: test with ambient in due time
            #ifdef MAQUIS_OPENMP
            #pragma omp parallel for collapse(2)
            #endif
		     	for (pos_t i = 0; i < lattice.size(); ++i){
				for (pos_t j = 0; j < lattice.size(); ++j){
					for (pos_t k = ((bra_neq_ket) ? 0 : std::min(i, j)); k < lattice.size(); ++k){
                                    if(i == j && i  == k)
                                         continue;
						for (pos_t l = 0; l < lattice.size(); ++l){

                                          if(i == l && k > j)
                                                continue;
                                          if(std::max(i,j)  < std::max(k,l)  && (j > l || l > k))
                                                continue;
                                          if(std::max(i,j)  > std::max(k,l)  && std::max(k,l) > j)
                                                continue;
                                          if(std::max(i,j) == std::max(k,l)  && i > l)
                                                continue;
                                          if(std::min(i,j) == std::min(k,l)  && std::max(i,j) <= std::max(k,l))
                                                continue;
                                          if((i == j && l  == i) || (i == k && l == i) || (j == k && (l == j || l > i)))
                                                continue;

							pos_t idx[] = { i,k,l,j };
							pos_t inv_count=0, n=4;
        					      for(pos_t c1 = 0; c1 < n - 1; c1++)
            					for(pos_t c2 = c1+1; c2 < n; c2++)
                					if(idx[c1] > idx[c2]) inv_count++;

							int phase = 1;
							if (inv_count % 2)
								phase = -1;

							std::vector<pos_op_t> op_string;
							op_string.push_back( std::make_pair(i, ops[0][0].first[lattice.get_prop<int>("type",i)]));
							op_string.push_back( std::make_pair(k, ops[0][1].first[lattice.get_prop<int>("type",k)]));
							op_string.push_back( std::make_pair(l, ops[0][2].first[lattice.get_prop<int>("type",l)]));
							op_string.push_back( std::make_pair(j, ops[0][3].first[lattice.get_prop<int>("type",j)]));

							////// CALL MPO MAKER /////////
							NTermsMPO<Matrix, SymmGroup> rdm_elem(lattice, identities, fillings, op_string, phase);
							MPO<Matrix, SymmGroup> mpo = rdm_elem.create_mpo();
							//////////////////////////////////

							//std::vector<typename MPS<Matrix, SymmGroup>::scalar_type> dct = multi_expval(bra_mps, ket_mps, mpo);
							typename MPS<Matrix, SymmGroup>::scalar_type dct = expval(bra_mps, ket_mps, mpo);
							
                                          //ecounter++;
                                           
							if(dct != 0.0) {
								//maquis::cout << std::fixed << std::setprecision(10) << i+1 << " " << j+1 << " " << k+1 << " " << l+1 << "\t" << dct << std::endl;
							
							
								std::vector<pos_t> label; label.push_back(i); label.push_back(j); label.push_back(k); label.push_back(l);
								std::vector<std::vector<pos_t> > num_labels;
								num_labels.push_back(label);
								std::vector<std::string> lbt = label_strings(lattice, num_labels);
								
								this->vector_results.push_back(dct);
								this->labels.push_back(lbt[0]);
							} 

						} // l loop
					} // k loop
				} // j loop
			} // i loop

			//maquis::cout << std::fixed << std::setprecision(10) << " number of elements measured: " << ecounter << std::endl;
        }
        
    private:
        Lattice lattice;
        positions_type positions_first;
        op_vec identities, fillings;
        std::vector<bond_element> ops;
        bool half_only, is_nn;
        std::string bra_ckp;
    };
}

#endif
