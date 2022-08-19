/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef OVERLAP_PROPAGATOR_H
#define OVERLAP_PROPAGATOR_H

#include <string>
#include <vector>
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/checks.h"

/**
 * @brief Class representing an MPS/MPS contraction. 
 * 
 * The tensor network to be contracted is now the following:
 * 
 *     o--o--o--o--o
 *     |  |  |  |  |
 *     x--x--x--x--x
 * 
 * which is as the contraction implemented by the [BoundaryPropagator]
 * object, without the MPO. The structure of the class reflects the 
 * philosophy of [BoundaryPropagator]. Two key differences are:
 * 
 * 1) that the method does not allow to access directly the "boundaries",
 *    but rather their partial contraction with an MPS, i.e.,
 *     
 *     o--o-- --o--o
 *     |  |  |  |  |
 *     x--x--x--x--x
 * 
 *    which is what is needed by the DMRG[ortho] method.
 * 
 * 2) multiple MPS tensors (x--x--x--x--x) are supported. Again, this is
 *    connected to the fact that in DMRG[ortho] one may need to orthogonalize
 *    wrt multiple MPSs.
 */
template<class Matrix, class SymmGroup, class Storage>
class OverlapPropagator {
public:
  // Types definition
  using MPSType = MPS<Matrix, SymmGroup>;
  using BlockMatrixType = block_matrix<typename storage::constrained<Matrix>::type, SymmGroup>;
  using Contraction = contraction::Engine<Matrix, typename storage::constrained<Matrix>::type, SymmGroup>;

  /** @brief Class constructor */
  OverlapPropagator(const MPSType& refMPS, const std::vector<std::string>& mpsFiles, BaseParameters & parms_,
                    int initSite=0)
    : refMPS_(refMPS), L_(refMPS_.length())
  {
    nOrthogonalMPSs_ = mpsFiles.size();
    for (int n = 0; n < nOrthogonalMPSs_; n++) {
      MPSType tmp;
      maquis::cout << "Orthogonal state " << n << " loaded from " << mpsFiles[n] << std::endl;
      maquis::checks::symmetry_check(parms_, mpsFiles[n]);
      maquis::checks::orbital_order_check(parms_, mpsFiles[n]);
      load(mpsFiles[n], tmp);
      orthoMPS_.emplace_back(std::move(tmp));
      maquis::checks::right_end_check(mpsFiles[n], orthoMPS_[n], refMPS_[L_-1].col_dim()[0].first);
      maquis::cout << "Right end: " << orthoMPS_[n][L_-1].col_dim() << std::endl;
    }
    initializeData(initSite);
  }

  /** @brief Class constructor from a vector of MPSs*/
  OverlapPropagator(const MPSType& refMPS, const std::vector<MPSType>& otherMPSs, int initSite=0)
    : refMPS_(refMPS), L_(refMPS_.length()), nOrthogonalMPSs_(otherMPSs.size()), orthoMPS_(otherMPSs)
  {
    initializeData(initSite);
  }

  /** @brief Contracts the boundary with an additional MPS */
  auto getOrthogonalVector(int iVector, int iSite) const {
    return contraction::site_ortho_boundaries(refMPS_[iSite], orthoMPS_[iVector][iSite],
                                              partialContractionLeft_[iVector][iSite], partialContractionRight_[iVector][iSite+1]);
  }

  /** @brief Propagation algorithm for the left overlap boundaries (see [BoundaryPropagator]) */
  inline void propagateLeftOverlapBoundaries(int siteInitial, int siteFinal) {
    if (siteInitial < siteFinal) {
      for (int iSite = siteInitial; iSite < siteFinal; iSite++) {
        for (int n = 0; n < nOrthogonalMPSs_; n++)
          partialContractionLeft_[n][iSite+1] = Contraction::overlap_left_step(refMPS_[iSite], orthoMPS_[n][iSite], partialContractionLeft_[n][iSite]);
      }
    }
  }

  /** @brief Propagation algorithm for the right boundaries (see [BoundaryPropagator]) */
  inline void propagateRightOverlapBoundaries(int siteInitial, int siteFinal) {
    if (siteInitial > siteFinal) {
      for (int iSite = siteInitial; iSite > siteFinal; iSite--) {
        for (int n = 0; n < nOrthogonalMPSs_; n++)
          partialContractionRight_[n][iSite] = Contraction::overlap_right_step(refMPS_[iSite], orthoMPS_[n][iSite], partialContractionRight_[n][iSite+1]);
      }
    }
  }

private:

  /** @brief Code that is shared by all constructors */
  void initializeData(int initSite) {
    // Prepares the data structure
    preparesDataStructure();
    // Fills the partial contractions
    generateLeftOverlapContractions(initSite);
    generateRightOverlapContractions(initSite);
  }

  /** @brief Prepares the data structure that accomondates the partial contraction */
  void preparesDataStructure() {
    partialContractionLeft_.resize(nOrthogonalMPSs_);
    partialContractionRight_.resize(nOrthogonalMPSs_);
    for (int n = 0; n < nOrthogonalMPSs_; n++) {
      partialContractionLeft_[n].resize(L_+1);
      partialContractionRight_[n].resize(L_+1);
      partialContractionLeft_[n][0] = refMPS_.left_boundary()[0];
      partialContractionRight_[n][L_] = refMPS_.right_boundary()[0];
    }
  }

  /** @brief Generates the left partial overlap contractions */
  void generateLeftOverlapContractions(int initSite) {
    propagateLeftOverlapBoundaries(0, initSite);
  }

  /** @brief Generates the right partial overlap contractions */
  void generateRightOverlapContractions(int initSite) {
    propagateRightOverlapBoundaries(L_-1, initSite);
  }

  // Class members
  const MPSType& refMPS_;
  std::vector<MPSType> orthoMPS_;
  int nOrthogonalMPSs_, L_;
  std::vector< std::vector<BlockMatrixType> > partialContractionLeft_, partialContractionRight_;
};

#endif // OVERLAP_PROPAGATOR_H