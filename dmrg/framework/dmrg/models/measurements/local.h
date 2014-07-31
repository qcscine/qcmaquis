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

#ifndef MEASUREMENTS_LOCAL_H
#define MEASUREMENTS_LOCAL_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/meas_prepare.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class local : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef generate_mpo::MPOMaker<Matrix, SymmGroup> generator;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
    public:
        
        local(std::string const& name_, const Lattice & lat,
              op_vec const & identities_, op_vec const & fillings_,
              std::vector<bond_element> const& terms)
        : base(name_)
        , lattice(lat)
        , identities(identities_)
        , fillings(fillings_)
        , is_bond(true)
        , mpo_terms(terms)
        {
            this->cast_to_real = all_true(mpo_terms.begin(), mpo_terms.end(), static_cast<bool (*)(bond_element const&)>(&is_hermitian_meas));
        }

        local(std::string const& name_, const Lattice & lat,
              op_vec const & identities_, op_vec const & fillings_,
              op_vec const& op)
        : base(name_)
        , lattice(lat)
        , identities(identities_)
        , fillings(fillings_)
        , is_bond(false)
        , site_term(op)
        , mpo_terms(std::vector<bond_element>(1, bond_element(1, std::make_pair(op, false)) ))
        {
            this->cast_to_real = is_hermitian_meas(site_term);
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results.clear();
            this->labels.clear();
            
            typedef typename SymmGroup::subcharge subcharge;
            if (!rmps || this->is_super_meas || is_bond)
                return evaluate_with_mpo(mps);
            
            int L = mps.size();
            this->vector_results.reserve(this->vector_results.size() + L);
            this->labels.reserve(this->labels.size() + L);
            
            for (typename Lattice::pos_t p = 0; p < L; ++p) {
                subcharge type = lattice.get_prop<subcharge>("type", p);
                typename Matrix::value_type res = 0.;
                bool evaluated = false;
                if (site_term[type].n_blocks() > 0) {
                    MPOTensor<Matrix, SymmGroup> temp;
                    temp.set(0, 0, site_term[type]);
                    
                    MPSTensor<Matrix, SymmGroup> vec2 = contraction::site_hamil2(mps[p], rmps.get().left(p), rmps.get().right(p), temp);
                    res += mps[p].scalar_overlap(vec2);
                    evaluated = true;
                }
                if (evaluated) {
                    this->vector_results.push_back(res);
                    this->labels.push_back( lattice.get_prop<std::string>("label", p) );
                }
            }
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new local(*this);
        }
        
        void evaluate_with_mpo(MPS<Matrix, SymmGroup> const& mps)
        {
            typedef typename SymmGroup::subcharge subcharge;
            typedef std::map<std::string, typename Matrix::value_type> result_type;
            result_type res;
            
            typename MPS<Matrix, SymmGroup>::scalar_type nn;
            if (this->is_super_meas)
                nn = dm_trace(mps, this->phys_psi);
            
            /// collect results from all mpo terms
            for (typename std::vector<bond_element>::const_iterator it = mpo_terms.begin(); it != mpo_terms.end(); ++it) {
                typedef std::map<std::string, MPO<Matrix, SymmGroup> > mpo_map;
                mpo_map mpos = meas_prepare::local(lattice, identities, fillings, *it);
                
                /// measure the value at each site / bond
                for (typename mpo_map::const_iterator mit = mpos.begin(); mit != mpos.end(); ++mit) {
                    typename result_type::iterator match = res.find(mit->first);
                    if (match == res.end())
                        boost::tie(match, boost::tuples::ignore) = res.insert( std::make_pair(mit->first, 0.) );
                    
                    if (!this->is_super_meas) {
                        match->second += expval(mps, mit->second);
                    } else {
                        MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mit->second, this->phys_psi);
                        // static_cast needed for icpc 12.x
                        typedef typename MPS<Matrix, SymmGroup>::scalar_type (*overlap_func)(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &);
                        typename MPS<Matrix, SymmGroup>::scalar_type val = ::overlap(super_mpo, mps);
                        match->second += val/nn;
                    }
                }
            }
            
            /// copy results to base
            this->vector_results.reserve(this->vector_results.size() + res.size());
            this->labels.reserve(this->labels.size() + res.size());
            for (typename result_type::const_iterator it = res.begin(); it != res.end(); ++it) {
                this->labels.push_back(it->first);
                this->vector_results.push_back(it->second);
            }
        }
        
    private:
        Lattice lattice;
        op_vec identities, fillings;
        bool is_bond;
        op_vec site_term;
        std::vector<bond_element> mpo_terms;
    };
    
}

#endif
