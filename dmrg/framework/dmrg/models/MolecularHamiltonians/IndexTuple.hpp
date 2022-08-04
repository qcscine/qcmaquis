/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2022 Laboratory for Physical Chemistry, ETH Zurich
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

#ifndef INDEX_TUPLE_H
#define INDEX_TUPLE_H

#include <array>
#include <iostream>
#include "dmrg/utils/align.h"

namespace chem {
namespace detail {

template<class SymmGroup, int NumberOfElements=4>
class IndexTuple
{
public:
  using ArrayType = std::array<int, NumberOfElements>;
  using InitializerType = std::initializer_list<int>;

  /** @brief Default class constructor - initialize every element to 0 */
  IndexTuple() {
    std::fill(tmpStorage.begin(), tmpStorage.end(), 0);
  }

  /** @brief Class constructor taking a series of indices as input */
  IndexTuple(InitializerType inputArray) {
    int idx = 0;
    for (const auto& iElement: inputArray) {
      tmpStorage[idx] = iElement;
      idx++;
    }
  }

  /** @brief Constructor from an array */
  IndexTuple(const ArrayType& inputArray) {
    tmpStorage = inputArray;
  }

  /** @brief Combination of two IndexTuple objects */
  template<int NumberOfInputElements>
  IndexTuple(const IndexTuple<SymmGroup, NumberOfInputElements>& firstTuple,
             const IndexTuple<SymmGroup, NumberOfInputElements>& secondTuple)
  {
    static_assert(2*NumberOfInputElements == NumberOfElements, "Combination of tuple does not match");
    for (int i = 0; i < NumberOfInputElements; i++) { 
      tmpStorage[i] = firstTuple[i]; 
      tmpStorage[i+4] = secondTuple[i]; 
    }
  }

  /** @brief Square bracket operator calls the underlying array function */
  const int& operator[](int i) const {
    return tmpStorage[i];
  }

  /** @brief Non-const overload */
  int& operator[](int i) {
    return tmpStorage[i];
  }

  /** In-place alignment */
  void align(bool isHermitian) {
    if (doAlignment)
      tmpStorage = maquis::detail::alignArray(tmpStorage, isHermitian);
  }

  /** @brief Non-const begin pointer */
  auto begin() { return tmpStorage.begin(); }

  /** @brief Const begin pointer */
  const auto begin() const { return tmpStorage.begin(); }

  /** @brief Non-const end pointer */
  auto end() { return tmpStorage.end(); }

  /** @brief Const end pointer */
  const auto end() const { return tmpStorage.end(); }

  /** @brief Getter for the underlying data */
  auto& data() { return tmpStorage; }
  const auto& data() const { return tmpStorage; }

  /** @brief Sign of the index (and of the underlying permutation */
  int sign() const {
    int inv_count = 0;
    for (int c1 = 0; c1 < NumberOfElements - 1; c1++)
        for (int c2 = c1+1; c2 < NumberOfElements; c2++)
            if (tmpStorage[c1] > tmpStorage[c2])
              inv_count++;
    return 1 - 2 * (inv_count % 2);
  }

private:
  ArrayType tmpStorage;
  static constexpr bool doAlignment = maquis::detail::AlignTraitClass<SymmGroup>::doRealign;
};

template<class SymmGroup, int NumberOfElements>
bool operator<(const IndexTuple<SymmGroup, NumberOfElements>& c1,
               const IndexTuple<SymmGroup, NumberOfElements>& c2)
{
    return c1.data() < c2.data();
}

} // namespace detail
} // namespace chem

#endif