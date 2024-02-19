/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014- Laboratory for Physical Chemistry, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#ifndef QC_CHEM_PARSE_INTEGRALS_H
#define QC_CHEM_PARSE_INTEGRALS_H

#include <algorithm>
#include <boost/archive/text_iarchive.hpp>
#include <boost/lambda/core.hpp>
#include <exception>
#include <filesystem>
#include <fstream>
#include <ios>
#include <iostream>
#include <istream>
#include <iterator>
#include <limits>
#include <memory>
#include <ostream>
#include <sstream>
#include <stdexcept>
#include <string>
#include <utility>
#include <vector>

#include "alps/numeric/matrix/matrix.hpp"
#include "dmrg/models/MolecularHamiltonians/IndexTuple.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/BaseParameters.h"
#include "integral_helper.h"
#include "integral_interface.h"
#include "parser_detail.h"
#include "utils/io.hpp"

namespace chem::detail {

/**
 * @brief Adds one element to the data structure used to store the Hamiltonian
 * terms.
 * @tparam T type associated with the sclar values.
 * @tparam SymmGroup Group describing the Hamiltonian symmetry.
 * @tparam IndexType Type used to store the indices.
 * @tparam HamiltonianType enum class indicating the type of Hamiltonian (-->
 * the number of indices in the FCIDUMP).
 * @param t pair <IndexType, T> associated with a single Hamiltonian term.
 * @param inv_order vector (for non-standard lattice mapping).
 * @param matrix_elements vector storing the Hamiltonian coefficients.
 * @param indices vector storing the indices corresponding to [matrix_elements].
 * @param do_align if true, align (i.e., put in canonical order) the index.
 * @param isHermitian (true if the Hamiltonian is Hermitian -- used to enforce
 * the orbital canonization)
 * @param cutoff positive number, cutoff used to neglect terms in the
 * Hamiltonian.
 */
template <
    class T, class SymmGroup, class IndexType, Hamiltonian HamiltonianType,
    HamiltonianTransformation Transcorrelated>
void updateIndices(
    const std::pair<IndexType, T>& t, const std::vector<int>& inv_order,
    std::vector<T>& matrix_elements, std::vector<IndexType>& indices,
    bool do_align, bool isHermitian, double cutoff
) {
  // Types declaration
  using TupleType = chem::detail::IndexTuple<
      SymmGroup, getIndexDim(HamiltonianType, Transcorrelated)>;
  // Functor class used to reorder the integral indices (used for custom
  // sorting)
  struct reorderer {
    int operator()(int p, const std::vector<int>& inv_order) {
      return p >= 0 ? inv_order[p] : p;
    }
  };
  // Actual function
  if (std::abs(t.second) > cutoff) {
    matrix_elements.push_back(t.second);
    TupleType tmp;
    int idx = 0;
    for (const auto& iElement : t.first) {
      tmp[idx] = reorderer()(iElement - 1, inv_order);
      idx++;
    }
    if (do_align) {
      tmp.align(isHermitian);
    }
    indices.push_back(tmp.data());
  }
}

struct FcidumpHeaderInfo {
  int norb;
  int nelec;
  int ms2;
};

static FcidumpHeaderInfo parse_header(std::istream& is) {
  int nelec;
  int norb;
  int ms2;

  std::string token;
  std::string line;
  while (std::getline(is, line)) {
    std::istringstream line_stream(line);
    while (std::getline(line_stream, token, '=')) {
      if (token.find("NORB") != std::string::npos) {
        line_stream >> norb;
      } else if (token.find("NELEC") != std::string::npos) {
        line_stream >> nelec;
      } else if (token.find("MS2") != std::string::npos) {
        line_stream >> ms2;
      }
    }
    if (token.find("&END") != std::string::npos) {
      break;
    }
  }
  return {norb, nelec, ms2};
}

/**
 * @brief Integral parser.
 * @tparam T type associated with the sclar values.
 * @tparam SymmGroup Group describing the Hamiltonian symmetry.
 * @param parms Parameter container
 * @param lat lattice object.
 * @param do_align if true, permutes the indices to have a common sorting.
 * @return std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T> >
 */
template <
    class T, class SymmGroup, Hamiltonian HamiltonianType,
    HamiltonianTransformation Transcorrelated =
        HamiltonianTransformation::Conventional>
inline std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T>>
parse_integrals(
    BaseParameters& parms, const Lattice& lat, bool do_align = true
) {
  // Types and variable definition
  using pos_t = Lattice::pos_t;
  using TupleType = chem::detail::IndexTuple<
      SymmGroup, getIndexDim(HamiltonianType, Transcorrelated)>;
  using IndexType = chem::index_type<HamiltonianType, Transcorrelated>;
  using IntegralTupleType = integral_tuple<T, HamiltonianType, Transcorrelated>;
  using IntegralMapType = integral_map<T, HamiltonianType, Transcorrelated>;
  static constexpr int numberOfIntegers =
      getIndexDim(HamiltonianType, Transcorrelated);
  static constexpr bool isHermitian =
      (Transcorrelated == HamiltonianTransformation::Conventional);
  std::string integralFileName =
      TranscorrelatedTraitClass<Transcorrelated>::getIntegralFileName();
  //
  std::vector<int> inv_order;
  std::vector<T> matrix_elements;
  alps::numeric::matrix<Lattice::pos_t> idx_;
  // Loads the ordering from input and determine inverse ordering.
  // Note that the inverse ordering is what is actually needed.
  // In fact, the input tells me which lattice size of the *new* order
  // corresponds to which lattice site of the *old* order. However, to convert
  // the FCIDUMP, one needs to know which site of the *new* order corresponds to
  // which site of the *old* one.
  std::vector<pos_t> order(lat.size());
  if (!parms.is_set("orbital_order")) {
    std::string s;
    for (pos_t p = 0; p < lat.size(); ++p) {
      order[p] = p + 1;
      s += (std::to_string(p + 1) + (p < (lat.size() - 1) ? "," : ""));
    }
    parms.set("orbital_order", s);
    // std::cout << "orbital order string " << s << std::endl;
  } else {
    order = parms["orbital_order"].as<std::vector<pos_t>>();
  }
  if (order.size() != lat.size()) {
    throw std::runtime_error(
        "orbital_order length is not the same as the number of orbitals\n"
    );
  }
  // The order starts with 1, so we remove 1 for coherence with C++ standards
  std::transform(
      order.begin(), order.end(), order.begin(), [](pos_t p) { return p - 1; }
  );
  inv_order.resize(order.size());
  for (int p = 0; p < order.size(); ++p) {
    inv_order[p] =
        std::distance(order.begin(), std::find(order.begin(), order.end(), p));
  }
  // == PARSING OF THE DATA ==
  std::vector<index_type<HamiltonianType, Transcorrelated>> indices;
  std::unique_ptr<std::istream> orb_string;
  // FCIDUMP integrals provided as a single string (undocumented, used only for
  // testing purposes) Note that, in this case, we don't expect any header.
  if (parms.is_set("integrals")) {
    std::string integrals = parms["integrals"];
    orb_string = std::make_unique<std::istringstream>(integrals);
  }
  // Integrals provided as a file
  else if (parms.is_set(integralFileName)) {
    maquis::cout << "Retrieving integrals from the file "
                 << parms[integralFileName] << std::endl;
    std::string integral_file = parms[integralFileName];
    if (!std::filesystem::exists(integral_file)) {
      throw std::runtime_error(
          "integral_file " + integral_file + " does not exist\n"
      );
    }
    orb_string = std::make_unique<std::ifstream>(integral_file.c_str());
    auto [norb, nelec, ms2] = parse_header(*orb_string);
    if (norb != lat.size()) {
      throw std::runtime_error(
          "The number of orbitals in the FCIDUMP (" + std::to_string(norb) +
          ") does not match the "
          "input file (" +
          std::to_string(lat.size()) + ")\n"
      );
    }
    if (nelec != parms["nelec"]) {
      std::cout << "!! WARNING: The number of electrons in the FCIDUMP ("
                << nelec << ") does not match the input file ("
                << parms["nelec"] << ") !!\n";
    }
  }
  // Integrals provided as a binary file
  else if (parms.is_set("integrals_binary")) {
    IntegralMapType ints;
    std::stringstream ss(parms["integrals_binary"].as<std::string>());
    boost::archive::text_iarchive ia{ss};
    ia >> ints;
    for (auto&& t : ints) {
      updateIndices<T, SymmGroup, IndexType, HamiltonianType, Transcorrelated>(
          t, inv_order, matrix_elements, indices, do_align, isHermitian,
          parms["integral_cutoff"]
      );
    }
  } else {
    throw std::runtime_error("Integrals are not defined in the input.");
  }

  // Read the FCIDUMP file/string and parse it. Only do it if the orb_string
  // pointer is not empty which is the case exactly when we want to parse the
  // FCIDUMP file (see above, i.e. when parms["integrals"] or
  // parms["integral_file"] is set. Otherwise, the pointer is empty, but
  // parms["integrals_binary"] is set and parsing is already completed, so the
  // below can be skipped. For this reason, there is a bit of code repetition
  // compared to above.

  if (orb_string) {
    T val;
    while (parser_detail::read_value<T>(*orb_string, val)) {
      IntegralTupleType t;
      t.second = val;
      // Parses the integral part.
      try {
        for (int iElement = 0; iElement < t.first.size(); iElement++) {
          *orb_string >> t.first[iElement];
        }
        //*(orb_string.get()) >> t.first[0] >> t.first[1] >> t.first[2] >>
        // t.first[3];
      } catch (std::exception& e) {
        std::cerr << e.what() << std::endl;
        throw std::runtime_error("error parsing integrals");
      }
      updateIndices<T, SymmGroup, IndexType, HamiltonianType, Transcorrelated>(
          t, inv_order, matrix_elements, indices, do_align, isHermitian,
          parms["integral_cutoff"]
      );
    }
  }

  // By now we should have parsed all the integrals, but we still have to
  // convert the indices to alps::numeric::matrix<Lattice::pos_t> Leon: I didn't
  // figure out how to safely add a row to alps::matrix using POD and not
  // iterators so I'm using a temporary object to read all the integrals and
  // then use resize on the alps::matrix once I know the temporary object's
  // size.

  idx_.resize(indices.size(), numberOfIntegers);

  // is better done with row iterators
  for (int i = 0; i < idx_.num_rows(); i++) {
    for (int j = 0; j < numberOfIntegers; j++) {
      idx_(i, j) = indices[i][j];
    }
  }

  // Integral dumping into HDF5 below MUST BE DISABLED if one builds
  // dmrg_multi_meas! Dumps the integrals into the result file for
  // reproducibility
  if (parms.is_set("donotsave") && parms["donotsave"] == 0 &&
      parms.is_set("resultfile")) {
    // dump indices but starting with 1 and with 0 as originally in the FCIDUMP
    std::vector<Lattice::pos_t> indices1;
    indices1.reserve(numberOfIntegers * indices.size());
    for (auto&& idx : indices) {
      for (auto&& i : idx) {
        indices1.push_back(i + 1);
      }
    }
    storage::archive ar(parms["resultfile"], "w");
    ar["/integrals/elements"] << matrix_elements;
    ar["/integrals/indices"] << indices1;
  }
// If in debug mode, checks that the indices are in the correct range.
#ifndef NDEBUG
  for (int m = 0; m < matrix_elements.size(); ++m)
    assert(
        *std::max_element(idx_.elements().first, idx_.elements().second) <=
        lat.size()
    );
#endif
  //
  return std::make_pair(idx_, matrix_elements);
}

}  // namespace chem::detail

#endif  // QC_CHEM_PARSE_INTEGRALS_H
