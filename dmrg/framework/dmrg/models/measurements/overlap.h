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

#ifndef MEASUREMENTS_OVERLAP_H
#define MEASUREMENTS_OVERLAP_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class overlap : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
    public:
        overlap(std::string const& name_,
                std::string const& ckp_)
        : base(name_)
        , bra_ckp(ckp_)
        { this->cast_to_real = false; }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            maquis::cout << "   overlap with " << bra_ckp << "." << std::endl;
            MPS<Matrix, SymmGroup> bra_mps;
            load(bra_ckp, bra_mps);
            
            if (bra_mps[bra_mps.length()-1].col_dim().sum_of_sizes() == 1)
                this->result = ::overlap(bra_mps, mps);
            else
                this->vector_results = multi_overlap(bra_mps, mps);
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new overlap(*this);
        }
        
    private:
        std::string bra_ckp;
    };
    
}

#endif
