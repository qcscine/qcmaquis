/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_NORMAL_ORDER_H
#define QC_NORMAL_ORDER_H

#include <array>
#include <vector>
#include <unordered_set>
#include <numeric>
#include <algorithm>

#include "dmrg/models/JordanWignerManager.h"

namespace {
    double calculatePermutationSign(const std::vector<int> &perm) {
        if (perm.size() == 1)
            return 1.;
        int transitionsCount = 0;
        for (int idx1 = 0; idx1 < perm.size() - 1; idx1++) {
            for (int idx2 = idx1 + 1; idx2 < perm.size(); idx2++) {
                if (perm[idx1] > perm[idx2])
                    transitionsCount += 1;
            }
        }
        return (transitionsCount % 2 == 0) ? 1. : -1.;
    }
}

/**
   * @brief Reorders the chain of operators and positions given in opVector and posVector
   *        to normal order, with respect to the hole states given in holeStates. The initial chain
   *        consists of only creation and annihilation operators. Returns permutation sign
 */
template<typename pos_t>
double applyNormalOrdering(std::vector<pos_t> &posVector, std::vector<OperatorType> &opVector,
                         const std::unordered_set<std::size_t> &holeStates) {
    const int K = posVector.size();

    // Create permutation vector, initialized 0,1,..,K-1
    std::vector<int> permutation(K);
    std::iota(permutation.begin(), permutation.end(), 0);

    std::stable_sort(permutation.begin(), permutation.end(), [&](int i, int j) {
        bool i_hole = holeStates.count(posVector[i]);
        bool j_hole = holeStates.count(posVector[j]);

        if (opVector[i] == OperatorType::CreateAlpha || opVector[i] == OperatorType::CreateBeta) {

            // Both creation
            if (opVector[j] == OperatorType::CreateAlpha || opVector[j] == OperatorType::CreateBeta) {
                if(i < j) {
                    // If i before j, ordering wrong iff i hole, j particle
                    return !(i_hole && !j_hole);
                } else {
                    // If j before i, order i before j iff i particle, j hole
                    return !i_hole && j_hole;
                }
            }

            // i is a Creation, j an Annihilation operator
            if (i < j) {
                // If i before j, ordering wrong iff both hole
                return !(i_hole && j_hole);
            } else {
                // If j < i, order i before j iff both particle
                return !i_hole && !j_hole;
            }
        } else {
            // Both are Annihilation
            if (opVector[j] == OperatorType::DestroyAlpha || opVector[j] == OperatorType::DestroyBeta) {
                if(i < j) {
                    // If i before j, ordering wrong iff j hole, i particle
                    return !(j_hole && !i_hole);
                } else {
                    // If j < i, order i before j iff j particle, i hole
                    return !j_hole && i_hole;
                }
            }

            // i is an Annihilation, j a Creation operator
            if (i < j) {
                // If i before j, ordering wrong iff both particle
                return !(!i_hole && !j_hole);
            } else {
                // If j < i, order i before j iff both hole
                return i_hole && j_hole;
            }
        }
    });

    double sign = calculatePermutationSign(permutation);

    std::vector<pos_t> tmpPosVector(K);
    std::vector<OperatorType> tmpOpVector(K);

    for (int i = 0; i < K; ++i) {
        tmpPosVector[i] = posVector[permutation[i]];
        tmpOpVector[i] = opVector[permutation[i]];
    }

    posVector = tmpPosVector;
    opVector = tmpOpVector;

    return sign;
}

#endif //QC_NORMAL_ORDER_H
