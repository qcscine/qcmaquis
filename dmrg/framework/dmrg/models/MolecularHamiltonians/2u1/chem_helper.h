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

#ifndef QC_CHEM_DETAIL_H
#define QC_CHEM_DETAIL_H

#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/MolecularHamiltonians/parse_integrals.h"
#include "dmrg/models/MolecularHamiltonians/2u1/term_maker.h"
#include "dmrg/models/OperatorHandlers/TagHandler.h"

namespace chem {
namespace detail {

template <
    typename Matrix, class SymmGroup, Hamiltonian HamiltonianType,
    HamiltonianTransformation Transcorrelated>
class ChemHelper {
 public:
  // Type definition
  // using value_type = typename Matrix::value_type;
  // Hardcoded to double because there is no reason to have complex
  // coefficients.
  using value_type = double;
  using term_descriptor = ::term_descriptor<value_type>;
  using tag_type = typename TagHandler<Matrix, SymmGroup>::tag_type;
  using pos_t = Lattice::pos_t;

  /** @brief Class constructor */
  ChemHelper(
      BaseParameters& parms, Lattice const& lat_,
      std::vector<tag_type> const& ident_, std::vector<tag_type> const& fill_,
      std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler_,
      bool doRealign = true
  )
      : lat(lat_), ident(ident_), fill(fill_), tag_handler(tag_handler_) {
    boost::tie(idx_, matrix_elements) = parse_integrals<
        value_type, SymmGroup, HamiltonianType, Transcorrelated>(
        parms, lat, doRealign
    );
    for (int m = 0; m < matrix_elements.size(); ++m) {
      IndexTuple<SymmGroup, numberOfIntegers> pos;
      std::copy(idx_.row(m).first, idx_.row(m).second, pos.begin());
      coefficients[pos] = matrix_elements[m];
    }
  }

  /** @brief Getter for a reference to the matrix elements */
  std::vector<value_type>& getMatrixElements() { return matrix_elements; }

  /** @brief Getter for the indices */
  int idx(int m, int pos) const { return idx_(m, pos); }

 private:
  static constexpr int numberOfIntegers =
      getIndexDim(HamiltonianType, Transcorrelated);
  const std::vector<tag_type>& ident;
  const std::vector<tag_type>& fill;
  std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
  const Lattice& lat;
  std::vector<value_type> matrix_elements;
  alps::numeric::matrix<Lattice::pos_t> idx_;
  std::map<IndexTuple<SymmGroup, numberOfIntegers>, value_type> coefficients;
};

}  // namespace detail
}  // namespace chem

#endif
