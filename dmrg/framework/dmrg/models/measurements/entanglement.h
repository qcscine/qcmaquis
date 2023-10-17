/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_ENTANGLEMENT_H
#define MEASUREMENTS_ENTANGLEMENT_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {

    template <class Matrix, class SymmGroup>
    class entropies : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
    public:
        entropies()
        : base("Entropy")
        {
            this->cast_to_real = true;
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results = calculate_bond_entropies(mps);
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new entropies(*this);
        }
    };
    
    template <class Matrix, class SymmGroup>
    class renyi_entropies : public measurement<Matrix, SymmGroup> {
        typedef  measurement<Matrix, SymmGroup> base;
    public:
        renyi_entropies()
        : base("Renyi2")
        {
            this->cast_to_real = true;
        }
        
        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            this->vector_results = calculate_bond_renyi_entropies(mps, 2);
        }
        
    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new renyi_entropies(*this);
        }
    };
    
}

#endif
