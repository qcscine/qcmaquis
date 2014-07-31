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

#ifndef MEASUREMENTS_AVERAGE_H
#define MEASUREMENTS_AVERAGE_H

#include "dmrg/models/measurement.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/super_mpo.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class average : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
        typedef generate_mpo::MPOMaker<Matrix, SymmGroup> generator;
        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
    public:
        
        /// Site term
        average(std::string const& name_,
                const Lattice & lattice,
                op_vec const & identities, op_vec const & fillings,
                op_vec const& ops)
        : base(name_)
        {
            generator mpom(lattice, identities, fillings);
            
            for (std::size_t p = 0; p < lattice.size(); ++p)
            {
                generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                term.operators.push_back( std::make_pair(p, ops[lattice.get_prop<int>("type", p)]) );
                mpom.add_term(term);
            }
            mpo = mpom.create_mpo();

            this->cast_to_real = is_hermitian_meas(ops);
        }
        
        /// Bond term
        average(std::string const& name_,
                const Lattice & lattice,
                op_vec const & identities, op_vec const & fillings,
                std::vector<bond_element> const& ops)
        : base(name_)
        {
            generator mpom(lattice, identities, fillings);
            
            for (std::size_t i = 0; i < ops.size(); ++i) {
                for (std::size_t p = 0; p < lattice.size(); ++p) {
                    generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                    term.operators.push_back( std::make_pair(p, ops[i][0].first[lattice.get_prop<int>("type", p)]) );
                    term.with_sign = ops[i][0].second;
                    std::vector<Lattice::pos_t> neighs = lattice.forward(p);
                    for (typename std::vector<Lattice::pos_t>::const_iterator hopto = neighs.begin();
                         hopto != neighs.end(); ++hopto)
                    {
                        generate_mpo::Operator_Term<Matrix, SymmGroup> term2(term);
                        term2.operators.push_back( std::make_pair(*hopto, ops[i][1].first[lattice.get_prop<int>("type", p)]) );
                        mpom.add_term(term2);
                    }
                }
            }
            mpo = mpom.create_mpo();
            
            /// TODO: this doesn't really work
            this->cast_to_real = all_true(ops.begin(), ops.end(), static_cast<bool (*)(bond_element const&)>(&is_hermitian_meas));
        }

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            if (!this->is_super_meas){
                this->result = expval(mps, mpo);
            } else {
                typename MPS<Matrix, SymmGroup>::scalar_type nn = dm_trace(mps, this->phys_psi);
                MPS<Matrix, SymmGroup> super_mpo = mpo_to_smps(mpo, this->phys_psi);
                // static_cast needed for icpc 12.x
#ifdef __INTEL_COMPILER
                typedef typename MPS<Matrix, SymmGroup>::scalar_type (*overlap_func)(MPS<Matrix, SymmGroup> const &, MPS<Matrix, SymmGroup> const &);
                this->result = static_cast<overlap_func>(&overlap)(super_mpo, mps) / nn;
#else
                this->result =::overlap<Matrix,SymmGroup>(super_mpo, mps) / nn;
#endif
            }
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new average(*this);
        }
        
    private:
        MPO<Matrix, SymmGroup> mpo;
    };
    
}

#endif
