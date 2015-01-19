/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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
    parse_integrals(BaseParameters & parms, Lattice const & lat)
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

        std::copy(order.begin(), order.end(), maquis::ostream_iterator<Lattice::pos_t>(std::cout, " "));
        maquis::cout << std::endl;
        std::copy(inv_order.begin(), inv_order.end(), maquis::ostream_iterator<Lattice::pos_t>(maquis::cout, " "));
        maquis::cout << std::endl;

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

        idx_.resize(raw.size()/5, 4);
        std::vector<double>::iterator it = raw.begin();
        int row = 0;
        while (it != raw.end()) {
            
            if (std::abs(*it) > parms["integral_cutoff"]){
                matrix_elements.push_back(*it++);
                std::vector<int> tmp;
                std::transform(it, it+4, std::back_inserter(tmp), boost::lambda::_1-1);

                IndexTuple aligned = align<SymmGroup>(reorderer()(tmp[0], inv_order), reorderer()(tmp[1], inv_order),
                                           reorderer()(tmp[2], inv_order), reorderer()(tmp[3], inv_order));
                idx_(row, 0) = aligned[0];
                idx_(row, 1) = aligned[1];
                idx_(row, 2) = aligned[2];
                idx_(row, 3) = aligned[3];
            }
            else { it++; }

            it += 4;
            row++;
        }

		/*   causing some troubles in the rel model!!!
        #ifndef NDEBUG
        for (std::size_t m = 0; m < matrix_elements.size(); ++m)
        {
            assert( *std::max_element(idx_.elements().first, idx_.elements().second) <= lat.size() );
        }
        #endif
		*/

        return std::make_pair(idx_, matrix_elements);
    }
}

#endif
