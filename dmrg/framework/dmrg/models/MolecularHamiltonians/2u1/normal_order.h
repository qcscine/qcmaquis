/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have receivedo a copy of the ALPS Application License along with
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
  if (perm.size() == 1) return 1.;
  int transitionsCount = 0;
  for (int idx1 = 0; idx1 < perm.size() - 1; idx1++) {
    for (int idx2 = idx1 + 1; idx2 < perm.size(); idx2++) {
      if (perm[idx1] > perm[idx2]) transitionsCount += 1;
    }
  }
  return (transitionsCount % 2 == 0) ? 1. : -1.;
}

std::string getIdentifier(int i, int j, int k, int l, int m, int n) {
  std::vector<std::pair<int, int>> indices = {{i, j}, {k, l}, {m, n}};

  // 2B is not hermitian, check if 2B
  if (!(i != -1 && j != -1 && k != -1 && l != -1 && m == -1 && n == -1)) {
    // Reorder indices in unique way (sorted)
    for (auto &v : indices)
      if (v.first > v.second) std::swap(v.first, v.second);
  }

  std::sort(indices.begin(), indices.end(), [](auto a, auto b) {
    return a.first < b.first || (a.first == b.first && a.second < b.second);
  });

  std::stringstream is;
  is << indices[0].first << " " << indices[0].second << " " << indices[1].first
     << " " << indices[1].second << " " << indices[2].first << " "
     << indices[2].second;

  return is.str();
}

std::string getIdentifier(int i, int j, int k, int l) {
  return getIdentifier(i, j, k, l, -1, -1);
}

std::string getIdentifier(int i, int j) {
  return getIdentifier(i, j, -1, -1, -1, -1);
}
}  // namespace

struct intPairHash {
  size_t operator()(std::pair<int, int> pair) const noexcept {
    return size_t(pair.first) << 32 | pair.second;
  }
};

struct intTupleHash {
  size_t operator()(std::tuple<int, int, int, int> t) const noexcept {
    size_t hash = std::hash<int>()(std::get<0>(t));
    hash ^= std::hash<int>()(std::get<1>(t)) + 0x9e3779b9 + (hash << 6) +
            (hash >> 2);
    hash ^= std::hash<int>()(std::get<2>(t)) + 0x9e3779b9 + (hash << 6) +
            (hash >> 2);
    hash ^= std::hash<int>()(std::get<3>(t)) + 0x9e3779b9 + (hash << 6) +
            (hash >> 2);
    return hash;
  }
};

/**
 * @brief Reorders the chain of operators and positions given in opVector and
 * posVector to normal order, with respect to the hole states given in
 * holeStates. The initial chain consists of only creation and annihilation
 * operators. Returns permutation sign
 */
template <typename pos_t>
double applyNormalOrdering(
    std::vector<pos_t> &posVector, std::vector<OperatorType> &opVector,
    const std::unordered_set<std::size_t> &holeStates
) {
  const int K = posVector.size();

  // Create permutation vector, initialized 0,1,..,K-1
  std::vector<int> permutation(K);
  std::iota(permutation.begin(), permutation.end(), 0);

  std::stable_sort(permutation.begin(), permutation.end(), [&](int i, int j) {
    bool i_hole = holeStates.count(posVector[i]);
    bool j_hole = holeStates.count(posVector[j]);

    if (opVector[i] == OperatorType::CreateAlpha ||
        opVector[i] == OperatorType::CreateBeta) {
      // Both creation
      if (opVector[j] == OperatorType::CreateAlpha ||
          opVector[j] == OperatorType::CreateBeta) {
        if (i < j) {
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
      if (opVector[j] == OperatorType::DestroyAlpha ||
          opVector[j] == OperatorType::DestroyBeta) {
        if (i < j) {
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

template <
    typename Matrix, typename SymmGroup, Hamiltonian HamiltonianType,
    HamiltonianTransformation Transcorrelated>
class NormalOrderingHelper {
 private:
  std::unordered_map<std::string, int> idxToMatrixElement;
  std::unordered_map<std::string, int> idxToMatrixElement_2B;

 public:
  NormalOrderingHelper(
      chem::detail::ChemHelper<
          Matrix, SymmGroup, HamiltonianType, Transcorrelated> &term_assistant
  ) {
    if (Transcorrelated != HamiltonianTransformation::Transcorrelated) return;

    for (int i_el = 0; i_el < term_assistant.getMatrixElements().size();
         ++i_el) {
      int i = term_assistant.idx(i_el, 0);
      int j = term_assistant.idx(i_el, 1);
      int k = term_assistant.idx(i_el, 2);
      int l = term_assistant.idx(i_el, 3);
      int m = term_assistant.idx(i_el, 4);
      int n = term_assistant.idx(i_el, 5);

      idxToMatrixElement[getIdentifier(i, j, k, l, m, n)] = i_el;
    }
  }

  /**
   * @brief Function that returns a map, corresponding to a chain a_p^\dag a_q
   * of one creation and one annihilation operator, where p and q are stored in
   * the pair, and gives the coefficient of the NO 3B contribution to this chain
   * of operators.
   */
  std::unordered_map<std::tuple<int, int, int, int>, double, intTupleHash>
  getNoTwoBodyCoefficients(
      chem::detail::ChemHelper<
          Matrix, SymmGroup, HamiltonianType, Transcorrelated> &term_assistant,
      std::size_t nelems, std::size_t L,
      const std::unordered_set<std::size_t> &hole_states
  ) {
    auto &matrix_elements = term_assistant.getMatrixElements();
    std::unordered_map<std::tuple<int, int, int, int>, double, intTupleHash>
        idxToCoefficient;
    std::vector<double> coeffs = {-1.0, 1.0 / 2.0, 1.0 / 2.0};
    for (auto i : hole_states) {
      for (int p = 0; p < L; ++p) {
        for (int q = 0; q < L; ++q) {
          // 3B contribution
          for (int r = 0; r < L; ++r) {
            for (int s = 0; s < L; ++s) {
              std::vector<std::string> matrixElementStrings = {
                  getIdentifier(i, i, p, s, q, r),
                  getIdentifier(i, s, p, i, q, r),
                  getIdentifier(i, r, p, s, q, i)};

              for (int n = 0; n < matrixElementStrings.size(); ++n) {
                if (!idxToMatrixElement.count(matrixElementStrings[n]))
                  continue;

                idxToCoefficient[{p, q, r, s}] +=
                    coeffs[n] *
                    matrix_elements
                        [idxToMatrixElement[matrixElementStrings[n]]];
              }
            }
          }
        }
      }
    }

    std::unordered_map<std::tuple<int, int, int, int>, double, intTupleHash>
        compressedMap;

    for (auto el : idxToCoefficient) {
      if (el.second != 0) compressedMap[el.first] = el.second;
    }

    return compressedMap;
  }

  /**
   * @brief Function that returns a map, corresponding to a chain a_p^\dag a_q
   * of one creation and one annihilation operator, where p and q are stored in
   * the pair, and gives the coefficient of the NO 3B contribution to this chain
   * of operators.
   */
  std::unordered_map<std::pair<int, int>, double, intPairHash>
  getNoOneBodyCoefficients(
      chem::detail::ChemHelper<
          Matrix, SymmGroup, HamiltonianType, Transcorrelated> &term_assistant,
      std::size_t nelems, std::size_t L,
      const std::unordered_set<std::size_t> &hole_states
  ) {
    auto &matrix_elements = term_assistant.getMatrixElements();

    std::unordered_map<std::pair<int, int>, double, intPairHash>
        idxToCoefficient;
    std::vector<double> coeffs_3B = {-1.0, 1.0, 2.0, -2.0};
    std::vector<double> coeffs_2B = {2.0, -1.0};

    for (int p = 0; p < L; ++p) {
      for (int q = 0; q < L; ++q) {
        for (auto i : hole_states) {
          // 2B contribution
          std::vector<std::string> matrixElementStrings_2B = {
              getIdentifier(i, i, p, q), getIdentifier(i, q, p, i)};
          for (int n = 0; n < matrixElementStrings_2B.size(); ++n) {
            if (!idxToMatrixElement.count(matrixElementStrings_2B[n])) continue;

            idxToCoefficient[{p, q}] +=
                coeffs_2B[n] *
                matrix_elements[idxToMatrixElement[matrixElementStrings_2B[n]]];
          }

          // 3B contribution
          for (auto j : hole_states) {
            std::vector<std::string> matrixElementStrings_3B = {
                getIdentifier(i, j, j, p, q, i),
                getIdentifier(i, j, j, i, p, q),
                getIdentifier(i, i, p, j, j, q),
                getIdentifier(i, i, j, j, p, q)};

            for (int n = 0; n < matrixElementStrings_3B.size(); ++n) {
              if (!idxToMatrixElement.count(matrixElementStrings_3B[n]))
                continue;

              idxToCoefficient[{p, q}] +=
                  coeffs_3B[n] *
                  matrix_elements
                      [idxToMatrixElement[matrixElementStrings_3B[n]]];
            }
          }
        }
      }
    }

    std::unordered_map<std::pair<int, int>, double, intPairHash> compressedMap;

    for (auto el : idxToCoefficient) {
      if (el.second != 0) compressedMap[el.first] = el.second;
    }

    return compressedMap;
  };

  /**
   * @brief Function that returns the zero-body contribution of the
   * normal-ordered 3B operator
   */
  double getNoZeroBodyContribution(
      chem::detail::ChemHelper<
          Matrix, SymmGroup, HamiltonianType, Transcorrelated> &term_assistant,
      std::size_t nelems, std::size_t L,
      const std::unordered_set<std::size_t> &hole_states
  ) {
    auto &matrix_elements = term_assistant.getMatrixElements();
    double contribution = 0;

    std::vector<double> coeffs_3B = {2.0, -2.0 / 3.0, -4.0 / 3.0};
    std::vector<double> coeffs_2B = {2.0, -1.0};
    std::vector<double> coeffs_1B = {2.0};

    for (auto i : hole_states) {
      // 1B contribution
      std::vector<std::string> matrixElementStrings_1B = {getIdentifier(i, i)};
      for (int n = 0; n < matrixElementStrings_1B.size(); ++n) {
        if (!idxToMatrixElement.count(matrixElementStrings_1B[n])) continue;

        contribution +=
            coeffs_1B[n] *
            matrix_elements[idxToMatrixElement[matrixElementStrings_1B[n]]];
      }

      for (auto j : hole_states) {
        // 2B contribution
        std::vector<std::string> matrixElementStrings_2B = {
            getIdentifier(i, i, j, j), getIdentifier(i, j, j, i)};

        for (int n = 0; n < matrixElementStrings_2B.size(); ++n) {
          if (!idxToMatrixElement.count(matrixElementStrings_2B[n])) continue;

          contribution +=
              coeffs_2B[n] *
              matrix_elements[idxToMatrixElement[matrixElementStrings_2B[n]]];
        }

        // 3B contribution
        for (auto k : hole_states) {
          std::vector<std::string> matrixElementStrings_3B = {
              getIdentifier(i, i, j, k, k, j), getIdentifier(i, k, j, i, k, j),
              getIdentifier(i, i, j, j, k, k)};

          for (int n = 0; n < matrixElementStrings_3B.size(); ++n) {
            if (!idxToMatrixElement.count(matrixElementStrings_3B[n])) continue;

            contribution +=
                coeffs_3B[n] *
                matrix_elements[idxToMatrixElement[matrixElementStrings_3B[n]]];
          }
        }
      }
    }
    return contribution;
  };
};

#endif  // QC_NORMAL_ORDER_H
