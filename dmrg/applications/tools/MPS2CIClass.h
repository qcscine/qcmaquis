/*****************************************************************************
 *
 * ALPS MPS DMRG Project
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

#include <cmath>
#include <iterator>
#include <iostream>
#include <string>
#include <sys/time.h>
#include <sys/stat.h>

#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/tokenizer.hpp>
#include <boost/tuple/tuple.hpp>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/matrix_types.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/chem/util.h"
#include "ci_encode.hpp"
#include "sampling.hpp"

/** @brief Calculator that manages the calculation of the overlap between an MPS and an ONV */
template<class Matrix, class SymmGroup>
class MPSCICalculator {
public:
  // Types definition
  typedef int pos_t;
  using MPSType = MPS<Matrix, SymmGroup>;
  using ChargeType = typename SymmGroup::subcharge;
  using DeterminantType = std::vector< std::vector<typename SymmGroup::charge> >;

  /** @brief class constructor */
  MPSCICalculator(std::string mpsFileName) {
    // Loads the MPS
    load(mpsFileName, mps);
    // Loads the corresponding parameters
    L = mps.length();
    auto Nup = mps[L-1].col_dim()[0].first[0];
    auto Ndown = mps[L-1].col_dim()[0].first[1];
    DmrgParameters parms;
    parms.set("site_types", chem::detail::infer_site_types(mps));
    // Extract physical basis for every site from MPS
    std::vector<ChargeType> irreps = parms["site_types"];
    phys_dims = chem::detail::make_2u1_site_basis<Matrix, SymmGroup>(L, Nup, Ndown, parms["site_types"]);
    for (pos_t q = 0; q < L; ++q)
        per_site.push_back(phys_dims[irreps[q]]);
  }

  /** @brief Calculates the overlap  with a bunch of input determinants */
  void calculateOverlap(std::string determinantName) const {
    // Loads the determinants
    auto determinants = parse_config<Matrix, SymmGroup>(std::string(), per_site);
    // printout the determinants
    for (pos_t q = 0; q < determinants.size(); ++q) {
      for (pos_t p = 0; p < L; ++p)
        std::cout << determinants[q][p];
      std::cout << std::endl;
    }
    // Set initial counter
    int i = 1;
    for (const auto& it: determinants) {
      maquis::cout << "CI coefficient of det " << i << " : " << extract_coefficient(mps, it) << std::endl;
      i++;
    }
    maquis::cout << std::endl;
  }

private:
  MPSType mps;
  std::vector<Index<SymmGroup> > per_site;
  std::vector<Index<SymmGroup> > phys_dims;
  int L;
};
