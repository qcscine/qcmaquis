/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
