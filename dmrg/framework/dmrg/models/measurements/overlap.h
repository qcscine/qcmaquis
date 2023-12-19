/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MEASUREMENTS_OVERLAP_H
#define MEASUREMENTS_OVERLAP_H

#include <utility>

#include "dmrg/models/measurement.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

namespace measurements {

/** @brief Measurement associated with the overlap with a fixed MPS */
template <class Matrix, class SymmGroup>
class overlap : public measurement<Matrix, SymmGroup> {
  using base =  measurement<Matrix, SymmGroup>;
  using MPSType = MPS<Matrix, SymmGroup>;
public:

  /** @brief Class constructor from a checkpoint file */
  overlap(const std::string& name, const std::string& checkPointFile) : base(name)
  { 
    this->cast_to_real = false;
    load(checkPointFile, mpsRef);
  }
        
  // Constructor
  overlap(const std::string& name, const MPSType& mpsInput) : base(name), mpsRef(mpsInput) 
  { 
    this->cast_to_real = false;
  }
  
  /** @brief Method to calculate the overlap */
  void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
  {
    this->result = ::overlap(mpsRef, mps);
  }
    
protected:
  measurement<Matrix, SymmGroup>* do_clone() const
  {
    return new overlap(*this);
  }
    
private:
  MPSType mpsRef;
};

} // namespace measurements

#endif
