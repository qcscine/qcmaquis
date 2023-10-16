/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_CUSTOM_H
#define MEASUREMENTS_CUSTOM_H

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {
    
    template <class Matrix, class SymmGroup>
    class custom : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
        typedef typename base::op_t op_t;
    public:
        custom(std::string const& name_,
               const Lattice & lat,
               std::vector<op_t> const & identities,
               std::vector<op_t> const & fillings,
               std::vector< std::vector< std::pair<int, op_t> > > const & ops)
        : base(name_)
        {
            generate_mpo::MPOMaker<Matrix, SymmGroup> mpom(lat, identities, fillings);
            
            for (int k = 0; k < ops.size(); ++k) {
                generate_mpo::OperatorTerm<Matrix, SymmGroup> term;
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
