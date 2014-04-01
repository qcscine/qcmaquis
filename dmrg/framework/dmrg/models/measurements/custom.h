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

#ifndef MEASUREMENTS_CUSTOM_H
#define MEASUREMENTS_CUSTOM_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class custom : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
    public:
        custom(std::string const& name_,
               const Lattice & lat,
               std::vector<block_matrix<Matrix, SymmGroup> > const & identities,
               std::vector<block_matrix<Matrix, SymmGroup> > const & fillings,
               std::vector< std::vector< std::pair<int, block_matrix<Matrix, SymmGroup> > > > const & ops)
        : base(name_)
        {
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
            
            for (int k = 0; k < ops.size(); ++k) {
                generate_mpo::Operator_Term<Matrix, SymmGroup> term;
                term.operators = ops[k];
                term.with_sign = true;
                mpom.add_term(term);
            }
            
            mpo = mpom.create_mpo();
            
            this->cast_to_real = false;
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->result = expval(mps, mpo);
        }
    
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new custom(*this);
        }
        
    private:
        MPO<Matrix, SymmGroup> mpo;
    };
    
}

#endif
