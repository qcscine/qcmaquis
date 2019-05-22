/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2014 Laboratory for Physical Chemistry, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_CHEM_PARSE_INTEGRALS_H
#define QC_CHEM_PARSE_INTEGRALS_H

#include "integral_interface.h"

namespace chem {
namespace detail {

    template <class T, class SymmGroup>
    inline // need inline as this will be compiled in multiple objects and cause linker errors otherwise
    std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<T> >
    parse_integrals(BaseParameters & parms, Lattice const & lat, bool do_align = true)
    {
        typedef Lattice::pos_t pos_t;

        std::vector<pos_t> inv_order;
        std::vector<T> matrix_elements;
        alps::numeric::matrix<Lattice::pos_t> idx_;

        struct reorderer
        {
            pos_t operator()(pos_t p, std::vector<pos_t> const & inv_order) {
                return p >= 0 ? inv_order[p] : p;
            }
        };

        // load ordering and determine inverse ordering
        std::vector<pos_t> order(lat.size());
        if (!parms.is_set("orbital_order"))
            for (pos_t p = 0; p < lat.size(); ++p)
                order[p] = p+1;
        else
            order = parms["orbital_order"].as<std::vector<pos_t> >();

        if (order.size() != lat.size())
            throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");

        std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);
        inv_order.resize(order.size());
        for (int p = 0; p < order.size(); ++p)
            inv_order[p] = std::distance(order.begin(), std::find(order.begin(), order.end(), p));

        // ********************************************************************
        // *** Parse orbital data *********************************************
        // ********************************************************************

        std::vector<index_type> indices;


        // Stream used to load FCIDUMP file/FCIDUMP string
        std::unique_ptr<std::istream> orb_string;

        if (parms.is_set("integrals")) // FCIDUMP integrals in a string
        {
            // if we provide parameters inline, we expect it to be in FCIDUMP format without the header
            std::string integrals = parms["integrals"];
            orb_string = std::unique_ptr<std::istringstream>(new std::istringstream(integrals));
        }
        else if (parms.is_set("integral_file")) // FCIDUMP file
        {
            std::string integral_file = parms["integral_file"];
            if (!boost::filesystem::exists(integral_file))
                throw std::runtime_error("integral_file " + integral_file + " does not exist\n");

            orb_string = std::unique_ptr<std::ifstream>(new std::ifstream(integral_file.c_str()));

            // ignore the FCIDUMP file header -- 1st four lines
            for (int i = 0; i < 4; ++i)
                orb_string.get()->ignore(std::numeric_limits<std::streamsize>::max(),'\n');

        }
        else if (parms.is_set("integrals_binary")) // Serialized integral object
        {
            // parse serialized integrals

            integral_map<T> ints;
            std::stringstream ss(parms["integrals_binary"].as<std::string>());

            boost::archive::text_iarchive ia{ss};
            ia >> ints;

            for (auto&& t: ints)
            {
                if (std::abs(t.second) > parms["integral_cutoff"])
                {
                    matrix_elements.push_back(t.second);
                    if (do_align)
                    {
                        IndexTuple aligned = align<SymmGroup>(reorderer()(t.first[0]-1, inv_order), reorderer()(t.first[1]-1, inv_order),
                                                reorderer()(t.first[2]-1, inv_order), reorderer()(t.first[3]-1, inv_order));
                        indices.push_back({ aligned[0], aligned[1], aligned[2], aligned[3] });
                    }
                    else
                        indices.push_back({ t.first[0]-1, t.first[1]-1, t.first[2]-1, t.first[3]-1 });
                }
            }
        }
        else
            throw std::runtime_error("Integrals are not defined in the input.");

        if (orb_string)
        // Read the FCIDUMP file/string and parse it. Only do it if the orb_string pointer is not empty
        // which is the case exactly when we want to parse the FCIDUMP file (see above, i.e. when
        // parms["integrals"] or parms["integral_file"] is set.
        // Otherwise, the pointer is empty,
        // but parms["integrals_binary"] is set and parsing is already completed, so the below can be skipped.
        {
            while(*(orb_string.get()))
            {
                integral_tuple<T> t;

                try {
                    *(orb_string.get()) >> t.second >> t.first[0] >> t.first[1] >> t.first[2] >> t.first[3];
                }
                catch(std::exception & e)
                {
                    std::cerr << e.what() << std::endl;
                    throw std::runtime_error("error parsing integrals");
                }
                // ignore integrals that are below the cutoff threshold
                if (std::abs(t.second) > parms["integral_cutoff"])
                {
                    matrix_elements.push_back(t.second);
                    if (do_align)
                    {
                        IndexTuple aligned = align<SymmGroup>(reorderer()(t.first[0]-1, inv_order), reorderer()(t.first[1]-1, inv_order),
                                                reorderer()(t.first[2]-1, inv_order), reorderer()(t.first[3]-1, inv_order));
                        indices.push_back({ aligned[0], aligned[1], aligned[2], aligned[3] });
                    }
                    else
                        indices.push_back({ t.first[0]-1, t.first[1]-1, t.first[2]-1, t.first[3]-1 });

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
        if (parms.is_set("donotsave") && parms["donotsave"] == 0 && parms.is_set("resultfile"))
        {
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

        #ifndef NDEBUG
        for (std::size_t m = 0; m < matrix_elements.size(); ++m)
        {
            assert( *std::max_element(idx_.elements().first, idx_.elements().second) <= lat.size() );
        }
        #endif

        return std::make_pair(idx_, matrix_elements);
    }

    // Template specialization for complex numbers and U1DG symmetry!!
    template <>
    inline // need inline as this will be compiled in multiple objects and cause linker errors otherwise
    std::pair<alps::numeric::matrix<Lattice::pos_t>, std::vector<std::complex<double> > >
    parse_integrals <std::complex<double>,U1DG> (BaseParameters & parms, Lattice const & lat, bool do_align)
    {
        typedef Lattice::pos_t pos_t;
        typedef std::complex<double> T;
        typedef U1DG SymmGroup;

        std::vector<pos_t> inv_order;
        std::vector<T> matrix_elements;
        alps::numeric::matrix<Lattice::pos_t> idx_;

        struct reorderer
        {
            pos_t operator()(pos_t p, std::vector<pos_t> const & inv_order) {
                return p >= 0 ? inv_order[p] : p;
            }
        };

        // load ordering and determine inverse ordering
        std::vector<pos_t> order(lat.size());
        if (!parms.is_set("orbital_order"))
            for (pos_t p = 0; p < lat.size(); ++p)
                order[p] = p+1;
        else
            order = parms["orbital_order"].as<std::vector<pos_t> >();

        if (order.size() != lat.size())
            throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");

        std::transform(order.begin(), order.end(), order.begin(), boost::lambda::_1-1);
        inv_order.resize(order.size());
        for (int p = 0; p < order.size(); ++p)
            inv_order[p] = std::distance(order.begin(), std::find(order.begin(), order.end(), p));

        // ********************************************************************
        // *** Parse orbital data *********************************************
        // ********************************************************************

        std::string integral_file = parms["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");

        std::ifstream orb_file;
        orb_file.open(integral_file.c_str());
        for (int i = 0; i < 4; ++i)
            orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');

        std::vector<double> raw;
        std::copy(std::istream_iterator<double>(orb_file), std::istream_iterator<double>(),
                    std::back_inserter(raw));

        idx_.resize(raw.size()/6, 4);
        std::vector<double>::iterator it = raw.begin();
        int row = 0;
        while (it != raw.end()) {

            // create complex number
            double re = *it++;
            double im = *it++;
            T integral_value(re, im);

            //DEBUG
            //maquis::cout << integral_value.real() << " " << integral_value.imag() << std::endl;

            if (std::abs(integral_value) > parms["integral_cutoff"]){
                matrix_elements.push_back(integral_value);
                std::vector<int> tmp;
                std::transform(it, it+4, std::back_inserter(tmp), boost::lambda::_1-1);

                if (do_align)
                {
                    IndexTuple aligned = align<SymmGroup>(reorderer()(tmp[0], inv_order), reorderer()(tmp[1], inv_order),
                                               reorderer()(tmp[2], inv_order), reorderer()(tmp[3], inv_order));
                    idx_(row, 0) = aligned[0];
                    idx_(row, 1) = aligned[1];
                    idx_(row, 2) = aligned[2];
                    idx_(row, 3) = aligned[3];
                }
                else
                {
                    idx_(row, 0) = tmp[0];
                    idx_(row, 1) = tmp[1];
                    idx_(row, 2) = tmp[2];
                    idx_(row, 3) = tmp[3];
                }
            }
            else { idx_.remove_rows(row--); }

            it += 4;
            row++;
        }

        // dump the integrals into the result file for reproducibility
        if (parms["donotsave"] == 0)
        {
            std::vector<T> m_;
            std::vector<Lattice::pos_t> i_;

            it = raw.begin();
            while (it != raw.end()) {
                double re = *it++;
                double im = *it++;
                T integral_value(re, im);

                m_.push_back(integral_value);
                std::copy(it, it+4, std::back_inserter(i_));
                it += 4;
            }

            storage::archive ar(parms["resultfile"], "w");
            ar["/integrals/elements"] << m_;
            ar["/integrals/indices"] << i_;
        }

        #ifndef NDEBUG
        for (std::size_t m = 0; m < matrix_elements.size(); ++m)
        {
            assert( *std::max_element(idx_.elements().first, idx_.elements().second) <= lat.size() );
        }
        #endif

        return std::make_pair(idx_, matrix_elements);
    }
}
}

#endif
