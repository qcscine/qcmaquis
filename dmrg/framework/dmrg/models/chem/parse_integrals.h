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

namespace chem_detail {

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
        if (!parms.is_set("orbital_order")){
            std::string s;
            for (pos_t p = 0; p < lat.size(); ++p){
                order[p] = p+1;
                s += (std::to_string(p+1)+ (p < (lat.size()-1) ? "," : ""));
            }
            //std::cout << "orbital order string " << s << std::endl;
            parms.set("orbital_order", s);
        }
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

        // Distinguish between integral files and string buffers

        std::unique_ptr<std::istream> orb_string;

        if (parms.is_set("integrals"))
        {
            // if we provide parameters inline, we expect it to be in FCIDUMP format without the header
            std::string integrals = parms["integrals"];
            orb_string = std::unique_ptr<std::istringstream>(new std::istringstream(integrals));
        }
        else if (parms.is_set("integral_file"))
        {
            std::string integral_file = parms["integral_file"];
            if (!boost::filesystem::exists(integral_file))
                throw std::runtime_error("integral_file " + integral_file + " does not exist\n");

            orb_string = std::unique_ptr<std::ifstream>(new std::ifstream(integral_file.c_str()));

            // ignore the FCIDUMP header -- 1st four lines
            for (int i = 0; i < 4; ++i)
                orb_string.get()->ignore(std::numeric_limits<std::streamsize>::max(),'\n');

        }
        else
            throw std::runtime_error("Integrals are not defined in the input.");

        std::vector<double> raw;
        std::copy(std::istream_iterator<double>(*(orb_string.get())), std::istream_iterator<double>(),
                    std::back_inserter(raw));

        idx_.resize(raw.size()/5, 4);
        std::vector<double>::iterator it = raw.begin();
        int row = 0;
        while (it != raw.end()) {

            if (std::abs(*it) > parms["integral_cutoff"]){
                matrix_elements.push_back(*it++);
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
            else { it++; idx_.remove_rows(row--); }

            it += 4;
            row++;
        }

        // dump the integrals into the result file for reproducibility
        if (parms.is_set("resultfile") && parms.is_set("donotsave") && parms["donotsave"] == 0)
        {
            std::vector<double> m_;
            std::vector<Lattice::pos_t> i_;

            it = raw.begin();
            while (it != raw.end()) {
                m_.push_back(*it++);
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

    // Template specialization for complex numbersi and U1DG symmetry!!
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
        if (!parms.is_set("orbital_order")){
            std::string s;
            for (pos_t p = 0; p < lat.size(); ++p){
                order[p] = p+1;
                s += (std::to_string(p+1)+ (p < (lat.size()-1) ? "," : ""));
            }
            parms.set("orbital_order", s);
            //std::cout << "orbital order string " << s << std::endl;
        }
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

#endif
