/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
 *               2017-2017 by Alberto Baiardi <abaiardi@phys.ethz.ch>
 *
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

namespace VibrationalIntegralParser {

/**
 * @brief Integral parser for potential operators encoded in the n-mode representation
 * 
 * Reads the FF from an input file and creates the Hamiltonian for the n-mode representation.
 * The potential is given in the following form:
 * i-n i-m           float_1   --> 1-body term (the first index is the mode, the second the basis set)
 * i-n i-m j-p j-q   float_2   --> 2-body term (again, the indexes of the mode are the same.
 * 
 * @param parms 
 * @param lat lattice
 * @param line_string string to be parsed 
 * @return std::pair< std::vector< std::size_t > , T > integer representation of the parsed string
 */
template <class T >
inline // need inline as this will be compiled in multiple objects and cause linker errors otherwise
std::pair< std::vector< std::size_t > , T >
NModeIntegralParser(BaseParameters & parms, Lattice const & lat, std::string line_string)
{
    using pos_t = Lattice::pos_t;
    std::pair< std::vector< post_t > , T >  result ;
    int max_nbody = parms["nmode_max_coupling"] ;
    // -- Main data parsing --
    std::vector< std::size_t > size_vec ;
    std::vector< std::string > line_splitted ;
    // Trim leading and final spaces in the string.
    line_string.erase(line_string.begin(), std::find_if(line_string.begin(), line_string.end(),
                                                        [&](int ch) { return !std::isspace(ch); } ) ) ;
    line_string.erase(std::find_if(line_string.rbegin(), line_string.rend(),
                                   [&](int ch) { return !std::isspace(ch); } ).base(), line_string.end()
    // Split the string into i-n, i-m, float_1 or i-n, i-m, j-p, j-q, float_2 etc etc
    boost::split(line_splitted, line_string, boost::is_any_of(" "), boost::token_compress_on) ;
    T scalar_res = atof(line_splitted[line_splitted.size()-1].c_str()) ;
    //TODO ALB Here can be greatly improved
    if (line_splitted.size() == 3) {
        if (max_nbody >= 1) {
            // One-body term
            std::vector<std::string> term1, term2;
            boost::split(term1, line_splitted[0], boost::is_any_of("-"));
            boost::split(term2, line_splitted[1], boost::is_any_of("-"));
            // Checks data consistency
            assert (term1.size() == 2 && term2.size() == 2);
            size_vec.push_back(std::stoul(term1[0]));
            size_vec.push_back(std::stoul(term1[1]));
            size_vec.push_back(std::stoul(term2[0]));
            size_vec.push_back(std::stoul(term2[1]));
            assert (size_vec[0] == size_vec[2]);
            result = std::make_pair(size_vec, scalar_res) ;
        }
    } else if (line_splitted.size() == 5) {
        if (max_nbody >= 2) {
            // Two-body term
            std::vector<std::string> term1, term2, term3, term4;
            boost::split(term1, line_splitted[0], boost::is_any_of("-"));
            boost::split(term2, line_splitted[1], boost::is_any_of("-"));
            boost::split(term3, line_splitted[2], boost::is_any_of("-"));
            boost::split(term4, line_splitted[3], boost::is_any_of("-"));
            // Checks data consistency
            assert (term1.size() == 2 && term2.size() == 2 && term3.size() == 2 && term4.size() == 2);
            size_vec.push_back(std::stoul(term1[0]));
            size_vec.push_back(std::stoul(term1[1]));
            size_vec.push_back(std::stoul(term2[0]));
            size_vec.push_back(std::stoul(term2[1]));
            size_vec.push_back(std::stoul(term3[0]));
            size_vec.push_back(std::stoul(term3[1]));
            size_vec.push_back(std::stoul(term4[0]));
            size_vec.push_back(std::stoul(term4[1]));
            assert (size_vec[0] == size_vec[2] && size_vec[4] == size_vec[6]);
            result = std::make_pair(size_vec, scalar_res) ;
        }
    } else if (line_splitted.size() == 7) {
        if (max_nbody >= 3) {
            // Three-body term
            std::vector<std::string> term1, term2, term3, term4, term5, term6;
            boost::split(term1, line_splitted[0], boost::is_any_of("-"));
            boost::split(term2, line_splitted[1], boost::is_any_of("-"));
            boost::split(term3, line_splitted[2], boost::is_any_of("-"));
            boost::split(term4, line_splitted[3], boost::is_any_of("-"));
            boost::split(term5, line_splitted[4], boost::is_any_of("-"));
            boost::split(term6, line_splitted[5], boost::is_any_of("-"));
            // Checks data consistency
            assert (term1.size() == 2 && term2.size() == 2 && term3.size() == 2 && term4.size() == 2 &&
                    term5.size() == 2 && term6.size() == 2);
            size_vec.push_back(std::stoul(term1[0]));
            size_vec.push_back(std::stoul(term1[1]));
            size_vec.push_back(std::stoul(term2[0]));
            size_vec.push_back(std::stoul(term2[1]));
            size_vec.push_back(std::stoul(term3[0]));
            size_vec.push_back(std::stoul(term3[1]));
            size_vec.push_back(std::stoul(term4[0]));
            size_vec.push_back(std::stoul(term4[1]));
            size_vec.push_back(std::stoul(term5[0]));
            size_vec.push_back(std::stoul(term5[1]));
            size_vec.push_back(std::stoul(term6[0]));
            size_vec.push_back(std::stoul(term6[1]));
            assert (size_vec[0] == size_vec[2] && size_vec[4] == size_vec[6] && size_vec[8] == size_vec[10]);
            result = std::make_pair(size_vec, scalar_res) ;
        }
    } else {
        throw std::runtime_error("Operators with more than three-body terms NYI") ;
    }
    return result ;
}

} // namespace VibrationalIntegralParser

#endif // DMRG_VIBRATIONAL

#endif
