/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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
#include "dmrg/SweepBasedAlgorithms/SweepOptimizationTypeTrait.h"
#include "dmrg/SweepBasedAlgorithms/SweepMPSContainer.h"

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

  /** @brief Class constructor from a single MPS */
  OverlapPropagator(const MPSType& refMPS, const MPSType& otherMPS, int initSite=0)
    : refMPS_(refMPS), L_(refMPS_.length()), nOrthogonalMPSs_(1)
  {
    orthoMPS_.push_back(otherMPS);
    initializeData(initSite);
  }

  /** @brief Class constructor from a vector of MPS */
  OverlapPropagator(const MPSType& refMPS, const std::vector<MPSType>& otherMPSs, int initSite=0)
    : refMPS_(refMPS), L_(refMPS_.length()), nOrthogonalMPSs_(otherMPSs.size()), orthoMPS_(otherMPSs)
  {
    initializeData(initSite);
  }

  /** @brief Contracts the boundary with an additional MPS */
  template<SweepOptimizationType SweepType>
  auto getOrthogonalVector(int iVector, int siteLeft, int siteRight) const {
    using SweepMPSContainerType = SweepMPSContainer<Matrix, SymmGroup, SweepType>;
    auto mpsContainer = SweepMPSContainerType(refMPS_);
    auto orthoMPSContainer = SweepMPSContainerType(orthoMPS_[iVector]);
    return contraction::site_ortho_boundaries(mpsContainer.getMPSTensor(siteLeft),
                                              orthoMPSContainer.getMPSTensor(siteLeft),
                                              partialContractionLeft_[iVector][siteLeft], partialContractionRight_[iVector][siteRight]);
  }

  /** @brief Overload where the index of the state is not specified */
  template<SweepOptimizationType SweepType>
  auto getOrthogonalVector(int siteLeft, int siteRight) const {
    return this->template getOrthogonalVector<SweepType>(0, siteLeft, siteRight);
  }

  /** @brief Getter for the number of orthogonal states */
  int getNumberOfOverlapMPSs() const {
    return nOrthogonalMPSs_;
  }

  /** @brief Propagation algorithm for the left overlap boundaries (see [BoundaryPropagator]) */
  inline void updateLeftOverlapBoundaries(int iSite) {
    if (iSite > 0 && iSite <= L_) {
      for (int n = 0; n < nOrthogonalMPSs_; n++)
        partialContractionLeft_[n][iSite] = Contraction::overlap_left_step(refMPS_[iSite-1], orthoMPS_[n][iSite-1], partialContractionLeft_[n][iSite-1]);
    }
  }

  /** @brief Propagation algorithm for the right boundaries (see [BoundaryPropagator]) */
  inline void updateRightOverlapBoundaries(int iSite) {
    if (iSite >= 0 && iSite < L_) {
      for (int n = 0; n < nOrthogonalMPSs_; n++)
        partialContractionRight_[n][iSite] = Contraction::overlap_right_step(refMPS_[iSite], orthoMPS_[n][iSite], partialContractionRight_[n][iSite+1]);
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
    for (int iSite = 1; iSite < initSite; iSite++)
      updateLeftOverlapBoundaries(iSite);
  }

  /** @brief Generates the right partial overlap contractions */
  void generateRightOverlapContractions(int initSite) {
    for (int iSite = L_-1; iSite > initSite; iSite--)
      updateRightOverlapBoundaries(iSite);
  }

  // Class members
  const MPSType& refMPS_;
  std::vector<MPSType> orthoMPS_;
  int nOrthogonalMPSs_, L_;
  std::vector< std::vector<BlockMatrixType> > partialContractionLeft_, partialContractionRight_;
};

#endif // OVERLAP_PROPAGATOR_H
