/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
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

#include "parser_detail.h"
#include "integral_interface.h"
#include "dmrg/utils/align.h"

namespace chem {
namespace detail {

/**
 * @brief Integral parser.
 * 
 * @tparam T type associated with the sclar values.
 * @tparam SymmGroup Group describing the Hamiltonian symmetry.
 * @param parms Parameter container
 * @param lat lattice object.
 * @param do_align if true, permutes the indices to have a common sorting.
 * @param numberOfIndices number of indices expected in the FCIDUMP file
 * (4 for conventional calculations, 6 for transcorrelated).
 * @return std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T> > 
 */
template <class T, class SymmGroup>
inline
std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T> >
parse_integrals(BaseParameters& parms, const Lattice& lat, bool do_align=true, bool numberOfIndices=4)
{
    // Types and variable definition
    using pos_t = Lattice::pos_t;
    std::vector<pos_t> inv_order;
    std::vector<T> matrix_elements;
    alps::numeric::matrix<Lattice::pos_t> idx_;

    // Functor class used to reorder the integral indices (used for custom sorting)
    struct reorderer
    {
        pos_t operator()(pos_t p, std::vector<pos_t> const & inv_order) {
            return p >= 0 ? inv_order[p] : p;
        }
    };

    // Loads the ordering from input and determine inverse ordering.
    // Note that the inverse ordering is what is actually needed.
    // In fact, the input tells me which lattice size of the *new* order corresponds to which
    // lattice site of the *old* order. However, to convert the FCIDUMP, one needs to know
    // which site of the *new* order corresponds to which site of the *old* one.
    std::vector<pos_t> order(lat.size());
    if (!parms.is_set("orbital_order")) {
        std::string s;
        for (pos_t p = 0; p < lat.size(); ++p) {
            order[p] = p+1;
            s += (std::to_string(p+1)+ (p < (lat.size()-1) ? "," : ""));
        }
        parms.set("orbital_order", s);
        //std::cout << "orbital order string " << s << std::endl;
    }
    else {
        order = parms["orbital_order"].as<std::vector<pos_t> >();
    }
    if (order.size() != lat.size())
        throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");
    // The order starts with 1, so we remove 1 for coherence with C++ standards
    std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);
    inv_order.resize(order.size());
    for (int p = 0; p < order.size(); ++p)
        inv_order[p] = std::distance(order.begin(), std::find(order.begin(), order.end(), p));

    // == PARSING OF THE DATA ==
    std::vector<index_type<Hamiltonian::Electronic>> indices;
    std::unique_ptr<std::istream> orb_string;
    // FCIDUMP integrals provided as a single string (undocumented, used only for testing purposes)
    // Note that, in this case, we don't expect any header.
    if (parms.is_set("integrals"))
    {
        std::string integrals = parms["integrals"];
        orb_string = std::unique_ptr<std::istringstream>(new std::istringstream(integrals));
    }
    // Integrals provided as a file
    else if (parms.is_set("integral_file")) {
        std::string integral_file = parms["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
        orb_string = std::unique_ptr<std::ifstream>(new std::ifstream(integral_file.c_str()));
        // Ignore the FCIDUMP file header -- 1st four lines
        for (int i = 0; i < 4; ++i)
            orb_string.get()->ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
    // Integrals provided as a binary file
    else if (parms.is_set("integrals_binary"))
    {
        integral_map<T> ints;
        std::stringstream ss(parms["integrals_binary"].as<std::string>());
        boost::archive::text_iarchive ia{ss};
        ia >> ints;
        for (auto&& t: ints) {
            if (std::abs(t.second) > parms["integral_cutoff"]) {
                matrix_elements.push_back(t.second);
                if (do_align) {
                    //using AlignerClass = maquis::detail::AlignerClass<maquis::detail::AlignTraitClass<SymmGroup>::doRealign>;
                    IndexTuple aligned = align<SymmGroup>(reorderer()(t.first[0]-1, inv_order), reorderer()(t.first[1]-1, inv_order),
                                                          reorderer()(t.first[2]-1, inv_order), reorderer()(t.first[3]-1, inv_order));
                    indices.push_back({ aligned[0], aligned[1], aligned[2], aligned[3] });
                }
                else
                    indices.push_back({ t.first[0]-1, t.first[1]-1, t.first[2]-1, t.first[3]-1 });
            }
        }
    }
    else {
        throw std::runtime_error("Integrals are not defined in the input.");
    }

    // Read the FCIDUMP file/string and parse it. Only do it if the orb_string pointer is not empty
    // which is the case exactly when we want to parse the FCIDUMP file (see above, i.e. when
    // parms["integrals"] or parms["integral_file"] is set.
    // Otherwise, the pointer is empty,
    // but parms["integrals_binary"] is set and parsing is already completed, so the below can be skipped.
    // For this reason, there is a bit of code repetition compared to above.

    if (orb_string) 
    {
        T val;
        while(parser_detail::read_value<T>(*(orb_string.get()), val)) {
            integral_tuple<T> t;
            t.second = val;
            // Parses the integral part.
            try {
                *(orb_string.get()) >> t.first[0] >> t.first[1] >> t.first[2] >> t.first[3];
            } 
            catch(std::exception & e) {
                std::cerr << e.what() << std::endl;
                throw std::runtime_error("error parsing integrals");
            }
                
            // Ignore integrals that are below the cutoff threshold
            if (std::abs(t.second) > parms["integral_cutoff"])
            {
                matrix_elements.push_back(t.second);
                if (do_align) {
                    IndexTuple aligned = align<SymmGroup>(reorderer()(t.first[0]-1, inv_order), reorderer()(t.first[1]-1, inv_order),
                                                          reorderer()(t.first[2]-1, inv_order), reorderer()(t.first[3]-1, inv_order));
                    indices.push_back({ aligned[0], aligned[1], aligned[2], aligned[3] });
                }
                else {
                    indices.push_back({ t.first[0]-1, t.first[1]-1, t.first[2]-1, t.first[3]-1 });
                }

            }
        }
    }
    
    // by now we should have parsed all the integrals, but we still have to convert the indices to alps::numeric::matrix<Lattice::pos_t>
    // Leon: I didn't figure out how to safely add a row to alps::matrix using POD and not iterators
    // so I'm using a temporary object to read all the integrals
    // and then use resize on the alps::matrix once I know the temporary object's size.
    
    idx_.resize(indices.size(), 4);

    // is better done with row iterators
    for (int i = 0; i < idx_.num_rows(); i++)
        for (int j = 0; j < 4; j++)
            idx_(i,j) = indices[i][j];

    // Integral dumping into HDF5 below MUST BE DISABLED
    // if one builds dmrg_multi_meas!
    // dump the integrals into the result file for reproducibility
    if (parms.is_set("donotsave") && parms["donotsave"] == 0 && parms.is_set("resultfile")) {
        // dump indices but starting with 1 and with 0 as originally in the FCIDUMP
        std::vector<Lattice::pos_t> indices1;
        indices1.reserve(4*indices.size());
        for (auto&& idx: indices)
            for (auto&& i: idx)
                indices1.push_back(i+1);
        storage::archive ar(parms["resultfile"], "w");
        ar["/integrals/elements"] << matrix_elements;
        ar["/integrals/indices"] << indices1;
    }

    // If in debug mode, checks that the indices are in the correct range.
    #ifndef NDEBUG
    for (int m = 0; m < matrix_elements.size(); ++m) {
        assert( *std::max_element(idx_.elements().first, idx_.elements().second) <= lat.size() );   
    }    
    #endif

    return std::make_pair(idx_, matrix_elements);
}

} // namespace detail
} // namespace chem

#endif
