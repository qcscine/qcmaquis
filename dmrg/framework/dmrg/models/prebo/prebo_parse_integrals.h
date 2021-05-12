//
// Created by robin on 11.05.21.
//

#ifndef MAQUIS_DMRG_PREBO_PARSE_INTEGRALS_H
#define MAQUIS_DMRG_PREBO_PARSE_INTEGRALS_H

#include "integral_interface.h"

namespace prebo {
namespace detail {

// now the integral parser can both real and complex integrals without specialization
template <class T, class SymmGroup>
inline
std::pair<std::vector<chem::index_type<chem::Hamiltonian::PreBO>>, std::vector<T> >
parse_integrals(BaseParameters & parms, Lattice const & lat)
{
    typedef Lattice::pos_t pos_t;

    std::vector<T> matrix_elements;

    // ********************************************************************
    // *** Parse orbital data *********************************************
    // ********************************************************************

    std::vector<chem::index_type<chem::Hamiltonian::PreBO>> indices;

    if (parms.is_set("integral_file")) {
        std::string integral_file = parms["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");

        std::ifstream orb_file;
        std::string line_string;
        //
        orb_file.open(integral_file.c_str());

        // Layout of the integral file
        // Particle type - Basis Function ... integral
        // One body terms
        // i-j     i-l     double                    --> Same particle type, same spin,
        //                                               but different basis
        //                                               functions are possible.
        // Two body terms
        // i-j     m-n     m-p     i-l     double    --> The ordering must follow
        //                                               the rule: a(i)+ a(m)+ a(m)- a(i)-
        //                                               The matrix element would be
        //                                               zero otherwise.
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
            // EMPTY LINE
            // if (line_splitted.size() == 0) continue;

            // Last value in string is assigned to the integral value
            integral = atof(line_splitted[line_splitted.size() - 1].c_str());

            chem::integral_tuple<T, chem::Hamiltonian::PreBO> t;
            t.second = integral;
            //
            // Nuclear repulsion
            // Set all indices to -1
            if (line_splitted.size() == 1) {
                std::fill(t.first.begin(), t.first.end(), -1);
            }
            else {
                // remove integral value from vector
                // now the vector contains all 2nd quantization operators.
                line_splitted.pop_back();

                assert(line_splitted.size() == 2 || line_splitted.size() == 4);
                std::vector<std::string> indices_str;
                // loop over all 2nd quant. operators in vector.
                for (const auto &sq_op_str : line_splitted) {
                    std::vector<std::string> temp;
                    boost::split(temp, sq_op_str, boost::is_any_of("-"));
                    indices_str.insert(indices_str.end(), std::make_move_iterator(temp.begin()),
                                       std::make_move_iterator(temp.end()));
                }
                assert(indices_str.size() == 4 || indices_str.size() == 8);
                for (auto i = 0; i < 8; ++i) {
                    if (i < indices_str.size()) {
                        t.first[i] = std::stoul(indices_str[i]);
                        assert(t.first[i] < lat.size());
                    } else
                        t.first[i] = -1;
                }
            }
            if (std::abs(t.second) > parms["integral_cutoff"]) {
                matrix_elements.push_back(t.second);
                indices.push_back(t.first);
            }
        }
    }
    else if (parms.is_set("integrals_binary")) { // Serialized integral object
        // parse serialized integrals
        chem::integral_map<T, chem::Hamiltonian::PreBO> ints;
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

    // Integral dumping into HDF5 below MUST BE DISABLED
    // if one builds dmrg_multi_meas!

    // dump the integrals into the result file for reproducibility
    if (parms.is_set("donotsave") && parms["donotsave"] == 0 && parms.is_set("resultfile"))
    {
        // dump indices but starting with 1 and with 0 as originally in the FCIDUMP
        std::vector<Lattice::pos_t> indices_vec;

        indices_vec.reserve(chem::getIndexDim(chem::Hamiltonian::PreBO)*indices.size());
        for (auto&& idx: indices)
            for (auto&& i: idx)
                indices_vec.push_back(i);

        storage::archive ar(parms["resultfile"], "w");
        ar["/integrals/elements"] << matrix_elements;
        ar["/integrals/indices"] << indices_vec;
    }

    return std::make_pair(indices, matrix_elements);
} // parse_integrals

} // namespace detail
} // namespace prebo


#endif //MAQUIS_DMRG_PREBO_PARSE_INTEGRALS_H
