/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
 *               2017-2017 by Alberto Baiardi <abaiardi@phys.ethz.ch>
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

#ifndef VIB_PARSE_INTEGRALS_H
#define VIB_PARSE_INTEGRALS_H

#ifdef DMRG_VIBRATIONAL

#include "integral_interface.h"
#include "VibrationalModelTraitClass.hpp"

namespace Vibrational {
namespace detail {

/**
 * @brief Integral parser for potential operators encoded in the n-mode representation
 * 
 * Reads the FF from an input file and creates the Hamiltonian for the n-mode representation.
 * The potential is given in the following form:
 * i-n i-        m           float_1   --> 1-body term (the first index is the mode, the second the basis set)
 * i-n i-m j-p j-q           float_2   --> 2-body term (again, the indexes of the mode are the same).
 * i-n i-m j-p j-q k-v k-w   float_3   --> 3-body term.
 * 
 * @param parms  Parameter container
 * @param lat DMRG lattice
 * @param line_string string to be parsed 
 * @return std::pair< std::vector< std::size_t > , T > integer representation of the parsed string
 */
template <class T>
inline std::pair<std::vector<chem::index_type<chem::Hamiltonian::VibrationalNMode>>, std::vector<T> >
NModeIntegralParser(BaseParameters & parms, Lattice const & lat)
{
    typedef Lattice::pos_t pos_t;
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
            T integral;
            // initialize splitted line
            std::vector<std::string> line_splitted;
            std::vector<std::size_t> size_vec;
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
            chem::integral_tuple<T, chem::Hamiltonian::VibrationalNMode> t;
            t.second = integral;
            // Remove integral value from vector
            line_splitted.pop_back();
            std::vector<std::string> indices_str;
            // loop over all 2nd quant. operators in vector.
            for (const auto &sq_op_str : line_splitted) {
                std::vector<std::string> temp;
                boost::split(temp, sq_op_str, boost::is_any_of("-"));
                indices_str.insert(indices_str.end(), std::make_move_iterator(temp.begin()),
                                    std::make_move_iterator(temp.end()));
            }
            assert(indices_str.size() == 4 || indices_str.size() == 8 || line_splitted.size() == 12);
            for (auto i = 0; i < 12; ++i) {
                if (i < indices_str.size()) {
                    t.first[i] = std::stoul(indices_str[i]);
                    assert(t.first[i] < lat.size());
                } else {
                    t.first[i] = -1;
                }
            }
            if (std::abs(t.second) > parms["integral_cutoff"]) {
                matrix_elements.push_back(t.second);
                indices.push_back(t.first);
            }
        }
    }
    // Serialized integral object
    else if (parms.is_set("integrals_binary")) {
        // parse serialized integrals
        chem::integral_map<T, chem::Hamiltonian::VibrationalNMode> ints;
        std::stringstream ss(parms["integrals_binary"].as<std::string>());
        boost::archive::text_iarchive ia{ss};
        ia >> ints;
        for (auto&& t: ints)
        {
            if (std::abs(t.second) > parms["integral_cutoff"])
            {
                matrix_elements.push_back(t.second);
                indices.push_back(t.first);
            }
        }
    }
    else
        throw std::runtime_error("Integrals are not defined in the input.");
    // dump the integrals into the result file for reproducibility
    if (parms.is_set("donotsave") && parms["donotsave"] == 0 && parms.is_set("resultfile"))
    {
        // dump indices but starting with 1 and with 0 as originally in the FCIDUMP
        std::vector<Lattice::pos_t> indices_vec;
        indices_vec.reserve(chem::getIndexDim(chem::Hamiltonian::VibrationalNMode)*indices.size());
        for (auto&& idx: indices)
            for (auto&& i: idx)
                indices_vec.push_back(i);
        storage::archive ar(parms["resultfile"], "w");
        ar["/integrals/elements"] << matrix_elements;
        ar["/integrals/indices"] << indices_vec;
    }
    return std::make_pair(indices, matrix_elements);
}

/**
 * @brief Integral file parser for the PES expressed as a Taylor expansion.
 * 
 * This routine expects to parse an integral file that is written in the following format:
 * 
 * coefficient    i   j   k   l   m   n 
 * 
 * Each line of the file is associated with a SQ operator expressed as:
 * 
 * coeff * o_i * o_j * o_k * o_l * o_m * o_n
 * 
 * where the operator o_i is:
 *  - (b_i^\dagger + b_i) if i > 0
 *  - i*(b_{-i}^\dagger - b_{-i}) if i < 0
 *  - the identity if i == 0
 * 
 * Note that the dimension of the array depends on the maximum allowed coupling degree, which is 
 * set in the TrivialGroup vibrational model trait class.
 * 
 * @tparam T scalar type associated with the Hamiltonian (real for most vibrational calculations)
 * @param parms parameter container
 * @param lat DMRG lattice object
 * @return std::vector< std::pair< std::array<int, maxCoupling>, T > > vector containing the coefficients
 */

template<class T>
inline std::vector< std::pair< std::array<int, VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings>, T > > 
    WatsonIntegralParser(BaseParameters & parms, Lattice const & lat) 
{
    // Types definition
    using pos_t = Lattice::pos_t;
    using KeyType = std::array<int, VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings>;
    using RetType = std::vector< std::pair< KeyType, T> > ;
    // Load ordering and determine inverse ordering
    std::vector<pos_t> inv_order;
    std::vector<pos_t> order(lat.size());
    if (!parms.is_set("sites_order"))
        for (pos_t p = 0; p < lat.size(); ++p)
            order[p] = p+1;
    else
        order = parms["sites_order"].as<std::vector<pos_t> >();
    if (order.size() != lat.size())
        throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");
    // Removes 1 (to fullfill the C++ convetion) and calculates the inverse map
    // (which is the one that is actually used in )
    std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);
    inv_order.resize(order.size());
    for (int p = 0; p < order.size(); ++p)
        inv_order[p] = std::distance(order.begin(), std::find(order.begin(), order.end(), p));
    // -- Parses orbital data --
    RetType ret;
    if (parms.is_set("integral_file")) {
        std::string integral_file = parms["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
        std::ifstream orb_file;
        orb_file.open(integral_file.c_str());
        std::vector<double> raw;
        std::copy(std::istream_iterator<double>(orb_file), std::istream_iterator<double>(),
                    std::back_inserter(raw));
        auto it = raw.begin();
        // Determines the maximum many-body coupling degree
        std::vector<bool> doCoupling(VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings, false);
        doCoupling[0] = true;
        int upperBound = (parms.is_set("watson_max_coupling")) ? parms["watson_max_coupling"] : 6;
        for (int iActive = 0; iActive < upperBound; iActive++)
            doCoupling[iActive] = true;
        // == Main loop ==
        int row = 0;
        while (it != raw.end()) {
            // Computes the coupling degree of the Hamiltonian term
            std::vector<int> tmp2(VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings, 0);
            std::vector<int>::iterator jnk_iter;
            std::copy(it+1, it+7, tmp2.begin());
            if (tmp2[2] == 0)
                tmp2.resize(2);
            else if (tmp2[3] == 0)
                tmp2.resize(3);
            else if (tmp2[4] == 0)
                tmp2.resize(4);
            else if (tmp2[5] == 0)
                tmp2.resize(5);
            std::sort(tmp2.begin(), tmp2.end());
            jnk_iter = std::unique(tmp2.begin(), tmp2.end());
            long coupl = std::distance(tmp2.begin(), jnk_iter);
            if (std::abs(*it) > parms["integral_cutoff"] && doCoupling[coupl-1] ) {
                T coefficient = *it++;
                KeyType tmp;
                for (int idx = 0; idx < VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings; idx++)
                    tmp[idx] = *(it+idx);
                for (int idx = 0; idx < VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings; idx++)
                    if (tmp[idx] > 0)
                        tmp[idx] = inv_order[tmp[idx]-1]+1;
                    else if (tmp[idx] < 0)
                        tmp[idx] = -inv_order[-tmp[idx]-1]-1;
                ret.push_back(std::make_pair(tmp, coefficient));
            }
            else {
                ++it;
            }
            it += VibrationalModelTraitClass<TrivialGroup>::maximumNumberOfCouplings;
            row++;
        }
    }
    else if (parms.is_set("integrals_binary")) {
        // parse serialized integrals
        chem::integral_map<T, chem::Hamiltonian::VibrationalCanonical> ints;
        std::stringstream ss(parms["integrals_binary"].as<std::string>());
        boost::archive::text_iarchive ia{ss};
        ia >> ints;
        for (auto&& t: ints)
            if (std::abs(t.second) > parms["integral_cutoff"])
                ret.push_back(std::make_pair(t.first, t.second));
    }
    return ret;
}

} // namespace detail
} // namespace vibrational

#endif // DMRG_VIBRATIONAL

#endif
