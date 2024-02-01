/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef VIBRONIC_PARSE_INTEGRALS_H
#define VIBRONIC_PARSE_INTEGRALS_H

//#ifdef DMRG_VIBRATIONAL

#include "integral_interface.h"
#include "dmrg/models/lattice/lattice.h"

namespace Vibrational {
namespace detail {

/**
 * @brief Class representing a line of an Vibronic Hamiltonian file
 * @tparam T Scalar type underyling the Hamiltonian (can be double or complex)
 */
template <class T>
struct term_input {
  // Types definition
  using index_type = int;

  /** @brief Class constructor */
  explicit term_input(std::string& input_line)
      : is_header(false), coefficient(0.), indexes(0), dimension(0) {
    std::vector<std::string> line_splitted;
    boost::trim_left(input_line);
    boost::trim_right(input_line);
    boost::split(
        line_splitted, input_line, boost::is_any_of(" "),
        boost::token_compress_on
    );
    // The header is identified by the string "el_st", followed by the
    // electronic states that are coupled.
    if (line_splitted[0] == "EL_ST") {
      is_header = true;
      indexes.push_back(std::stoi(line_splitted[1]));
      indexes.push_back(std::stoi(line_splitted[2]));
    } else {
      coefficient = atof(line_splitted[0].c_str());
      for (std::size_t idx = 1; idx < line_splitted.size(); idx++)
        indexes.push_back(std::stoi(line_splitted[idx]));
      indexes.erase(
          std::remove(indexes.begin(), indexes.end(), 0), indexes.end()
      );
      dimension = indexes.size();
    }
  };

  /** True if the line is a header */
  bool is_header;
  /** Indices associated with the Hamiltonian term */
  std::vector<index_type> indexes;
  /** Scalar coefficient */
  T coefficient;
  /** Dimension of the Hamiltonian term */
  int dimension;
};

/**
 * @brief Parser for a vibronic Hamiltonian
 * @param parms parameters container
 * @param lat lattice
 * @return std::pair<std::vector<chem::index_type<chem::Hamiltonian::Vibronic>>,
 * std::vector<T> > Pair where the first element identifies the term of the
 * Hamiltonian and the second one identifies the coefficient.
 */

template <class T>
inline std::pair<
    std::vector<chem::index_type<chem::Hamiltonian::Vibronic>>, std::vector<T>>
parseIntegralVibronic(BaseParameters& parms, const Lattice& lat) {
  // Types and variables definition
  using pos_t = Lattice::pos_t;
  std::vector<T> matrix_elements;
  std::vector<chem::index_type<chem::Hamiltonian::Vibronic>> indices;
  // Variable initialization
  int n_states = parms["vibronic_num_elestates"];
  int L_lattice = lat.size();
  // == Parses orbital data ==
  std::string integral_file = parms["integral_file"];
  if (!std::filesystem::exists(integral_file))
    throw std::runtime_error(
        "integral_file " + integral_file + " does not exist\n"
    );
  std::ifstream orb_file;
  orb_file.open(integral_file.c_str());
  // -- MAIN LOOP --
  int state_1 = 0, state_2 = 0;
  std::string tmp;
  bool is_init = true;
  while (std::getline(orb_file, tmp)) {
    chem::integral_tuple<T, chem::Hamiltonian::Vibronic> t;
    // Computes the coupling degree of the Hamiltonian term
    auto line_obj = term_input<T>(tmp);
    // Check if the first element is a header
    if (is_init) {
      is_init = false;
      if (!line_obj.is_header)
        throw std::runtime_error("First line of the VibHam file must be header"
        );
    }
    // Check if it is a header
    if (line_obj.is_header) {
      state_1 = line_obj.indexes[0];
      state_2 = line_obj.indexes[1];
      // Otherwise it is a term
    } else {
      t.first[0] = state_1;
      t.first[1] = state_2;
      t.first[2] = line_obj.indexes[0];
      t.first[3] = line_obj.indexes[1];
      t.second = line_obj.coefficient;
      matrix_elements.push_back(t.second);
      indices.push_back(t.first);
    }
  }
  orb_file.close();
  return std::make_pair(indices, matrix_elements);
}

/**
 * @brief Parser for the integral file associated with the excitonic HH
 * Hamiltonian.
 *
 * The integral file is expected to be given in the following format:
 *
 *  i   i  coeff --> harmonic potential term
 * ...
 * -i  -i  coeff --> harmonic kinetic term
 * ...
 *  i   0  coeff --> LVC coupling term
 *
 * @tparam T Scalar Type underlying the definition of the wave function.
 * @param parms Parameter container
 * @param lat DMRG lattice
 * @return std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T> >
 */
template <class T>
inline std::pair<
    std::vector<chem::index_type<chem::Hamiltonian::Excitonic>>, std::vector<T>>
parseIntegralExcitonic(BaseParameters& parms, const Lattice& lat) {
  // Types and variables definition
  using pos_t = Lattice::pos_t;
  std::vector<T> matrix_elements;
  std::vector<chem::index_type<chem::Hamiltonian::Excitonic>> indices;
  std::string integral_file = parms["integral_file"];
  if (!std::filesystem::exists(integral_file))
    throw std::runtime_error(
        "integral_file " + integral_file + " does not exist\n"
    );
  //
  std::ifstream orb_file;
  orb_file.open(integral_file.c_str());
  std::vector<double> raw;
  std::copy(
      std::istream_iterator<double>(orb_file), std::istream_iterator<double>(),
      std::back_inserter(raw)
  );
  auto it = raw.begin();
  while (it != raw.end()) {
    // Computes the coupling degree of the Hamiltonian term
    std::vector<int> tmp;
    chem::integral_tuple<T, chem::Hamiltonian::Excitonic> t;
    t.second = *it;
    it++;
    std::copy(it, it + 2, std::back_inserter(tmp));
    t.first[0] = tmp[0];
    t.first[1] = tmp[1];
    it += 2;
    matrix_elements.push_back(t.second);
    indices.push_back(t.first);
  }
  return std::make_pair(indices, matrix_elements);
}

// VAL(i): extended Excitonic Hamiltonian parser
template <class T>
inline std::pair<std::vector<std::vector<int>>, std::vector<T>>
parseIntegralExcitonicExtended(BaseParameters& parms, const Lattice& lat) {
  // Types and variables definition
  using pos_t = Lattice::pos_t;
  std::vector<T> matrix_elements;
  std::string integral_file = parms["integral_file"];
  if (!boost::filesystem::exists(integral_file))
    throw std::runtime_error(
        "integral_file " + integral_file + " does not exist\n"
    );
  std::ifstream orb_file;
  orb_file.open(integral_file.c_str());
  // -- MAIN LOOP --
  std::string tmp;
  std::vector<std::string> line_splitted;
  std::vector<int> indices_tmp;
  std::vector<std::vector<int>> indices;
  while (std::getline(orb_file, tmp)) {
    boost::trim_left(tmp);
    boost::trim_right(tmp);
    boost::split(
        line_splitted, tmp, boost::is_any_of(" "), boost::token_compress_on
    );  // be careful with tabs. use spaces, otherwise FCIDUMP file is not read
        // in properly
    double coefficient = atof(line_splitted[0].c_str());
    for (std::size_t idx = 1; idx < line_splitted.size(); idx++) {
      indices_tmp.push_back(std::stoi(line_splitted[idx]));
    }
    indices.push_back(indices_tmp);
    matrix_elements.push_back(coefficient);
    indices_tmp.clear();
  }
  return std::make_pair(indices, matrix_elements);
}


template<class T>
inline std::pair<std::vector<chem::index_type<chem::Hamiltonian::VibrationalNMode>>, std::vector<T> >
    parseIntegralExcitonicNmode(BaseParameters& parms, const Lattice& lat)
{
    typedef Lattice::pos_t pos_t;
    using InputType = double;
    std::vector<T> matrix_elements;
    std::vector<chem::index_type<chem::Hamiltonian::VibrationalNMode>> indices;
    // Data read from input file
    if (parms.is_set("integral_file")) {
        std::string integral_file = parms["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
        std::ifstream orb_file;
        std::string line_string;
        orb_file.open(integral_file.c_str());
        // File parsing
        while (getline(orb_file, line_string)) {
            if (line_string[0] == '#' || line_string == "")
                continue;
            // initialize integral value
            InputType integral;
            // initialize splitted line
            std::vector<std::string> line_splitted;
            // -- Main data parsing --
            // Trim leading and final spaces in the string.
            line_string.erase(line_string.begin(),
                              std::find_if(line_string.begin(), line_string.end(),
                                           [&](int ch) { return !std::isspace(ch); }));
            line_string.erase(
                    std::find_if(line_string.rbegin(), line_string.rend(),
                                 [&](int ch) { return !std::isspace(ch); }).base(),
                    line_string.end());
            // Split the string
            boost::split(line_splitted, line_string, boost::is_any_of(" "), boost::token_compress_on);
            // Last value in string is assigned to the integral value
            integral = atof(line_splitted[line_splitted.size() - 1].c_str());
            chem::integral_tuple<InputType, chem::Hamiltonian::VibrationalNMode> t;
            t.second = integral;
            // Remove integral value from vector
            line_splitted.pop_back();
            std::vector<std::string> indices_str;
            // loop over all 2nd quant. operators in vector.
            for (const auto &sq_op_str : line_splitted) {
                std::vector<std::string> temp;
                boost::split(temp, sq_op_str, boost::is_any_of("-\t ")); //splits string whenever a hyphen, tab or space is encountered 
                indices_str.insert(indices_str.end(), std::make_move_iterator(temp.begin()),
                                    std::make_move_iterator(temp.end()));
            }
            assert(indices_str.size() == 4+2 || indices_str.size() == 8+2 || indices_str.size() == 12+2); //+2 as there si a number denoting excited state and whether mode connects two monomers
            for (auto i = 0; i < 12; ++i) {
                if (i < indices_str.size()) {
                    t.first[i] = std::stoul(indices_str[i]);
                } else {
                    t.first[i] = -1;
                }
            }
            if (std::abs(t.second) > parms["integral_cutoff"]) {
                matrix_elements.push_back(static_cast<T>(t.second));
                indices.push_back(t.first);
            }
        }
        orb_file.close();
    }
    else
        throw std::runtime_error("Integrals are not defined in the input.");
    return std::make_pair(indices, matrix_elements);

}
//VAL(f)

}  // namespace detail
}  // namespace Vibrational

// #endif // DMRG_VIBRATIONAL

#endif  // VIB_PARSE_INTEGRALS_H
