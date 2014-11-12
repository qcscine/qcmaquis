/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_CHEM_HELPER_H
#define QC_CHEM_HELPER_H

namespace chem_detail {

    template <typename M, class S>
    class ChemHelper
    {
    public:
        typedef typename M::value_type value_type;
        typedef ::term_descriptor<value_type> term_descriptor;
        typedef typename TagHandler<M, S>::tag_type tag_type;
        typedef Lattice::pos_t pos_t;

        ChemHelper(BaseParameters & parms, Lattice const & lat,
                   tag_type ident_, tag_type fill_, boost::shared_ptr<TagHandler<M, S> > tag_handler_) 
            : ident(ident_), fill(fill_), tag_handler(tag_handler_)
        {
            this->parse_integrals(parms, lat);

            for (std::size_t m=0; m < matrix_elements.size(); ++m) {
                IndexTuple pos;
                std::copy(idx_.row(m).first, idx_.row(m).second, pos.begin());
                coefficients[pos] = matrix_elements[m];
            }
        }

        std::vector<value_type> & getMatrixElements() { return matrix_elements; }
        
        int idx(int m, int pos) const {
            //return idx_[m][pos];
            return idx_(m,pos);
        }

        void commit_terms(std::vector<term_descriptor> & tagterms) {
            for (typename std::map<IndexTuple, term_descriptor>::const_iterator it = two_terms.begin();
                    it != two_terms.end(); ++it)
                tagterms.push_back(it->second);

            for (typename std::map<TermTuple, term_descriptor>::const_iterator it = three_terms.begin();
                    it != three_terms.end(); ++it)
                tagterms.push_back(it->second);
        }

        void add_term(std::vector<term_descriptor> & tagterms,
                      value_type scale, int p1, int p2, tag_type op_1, tag_type op_2) {

            term_descriptor
            term = TermMaker<M, S>::two_term(false, ident, scale, p1, p2, op_1, op_2, tag_handler);
            IndexTuple id(p1, p2, op_1, op_2);
            if (two_terms.count(id) == 0) {
                two_terms[id] = term;
            }
            else 
                two_terms[id].coeff += term.coeff;
        }

        void add_term(std::vector<term_descriptor> & tagterms,
                      value_type scale, int s, int p1, int p2, tag_type op_i, tag_type op_k, tag_type op_l, tag_type op_j) {

            term_descriptor
            term = TermMaker<M, S>::three_term(ident, fill, scale, s, p1, p2, op_i, op_k, op_l, op_j, tag_handler);
            TermTuple id(IndexTuple(s,s,p1,p2),IndexTuple(op_i,op_k,op_l,op_j));
            if (three_terms.count(id) == 0) {
                three_terms[id] = term;
            }
            else
                three_terms[id].coeff += term.coeff;
    
        }

        void add_term(std::vector<term_descriptor> & tagterms,
                      int i, int k, int l, int j, tag_type op_i, tag_type op_k, tag_type op_l, tag_type op_j)
        {
            // Collapse terms with identical operators and different scales into one term
            if (op_i == op_k && op_j == op_l) {

                // if i>j, we switch l,j to get the related term
                // if j<i, we have to switch i,k, otherwise we get a forbidden permutation
                IndexTuple self(i,j,k,l), twin(i,l,k,j);
                if (i<j) twin = IndexTuple(k,j,i,l);

                if (self > twin) {
                
                    term_descriptor
                    term = TermMaker<M, S>::four_term(ident, fill, coefficients[align(i,j,k,l)], i,k,l,j,
                                                   op_i, op_k, op_l, op_j, tag_handler);

                    term.coeff += value_type(sign(twin)) * coefficients[align(twin)];

                    tagterms.push_back(term);
                }
                //else: we already have the term
            }
            else {
                tagterms.push_back( TermMaker<M, S>::four_term(ident, fill, coefficients[align(i,j,k,l)], i,k,l,j,
                                   op_i, op_k, op_l, op_j, tag_handler) );
            }
        }
    
    private:
        void parse_integrals(BaseParameters & parms, Lattice const & lat) {

            // load ordering and determine inverse ordering
            std::vector<pos_t> order(lat.size());
            if (!parms.is_set("orbital_order"))
                for (pos_t p = 0; p < lat.size(); ++p)
                    order[p] = p+1;
            else
                order = parms["orbital_order"].template as<std::vector<pos_t> >();

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

                    IndexTuple aligned = align(reorder(tmp[0]), reorder(tmp[1]), reorder(tmp[2]), reorder(tmp[3]));
                    idx_(row, 0) = aligned[0];
                    idx_(row, 1) = aligned[1];
                    idx_(row, 2) = aligned[2];
                    idx_(row, 3) = aligned[3];
                }
                else { it++; }

                it += 4;
                row++;
            }

            #ifndef NDEBUG
            for (std::size_t m = 0; m < matrix_elements.size(); ++m)
            {
                assert( *std::max_element(idx_.elements().first, idx_.elements().second) <= lat.size() );
            }
            #endif
            
        }

        int reorder(int p) {
            return p >= 0 ? inv_order[p] : p;
        }

        tag_type ident, fill;
        boost::shared_ptr<TagHandler<M, S> > tag_handler;

        std::vector<value_type> matrix_elements;
        alps::numeric::matrix<Lattice::pos_t> idx_;
        std::vector<int> order;
        std::vector<int> inv_order;

        std::map<IndexTuple, value_type> coefficients;

        std::map<TermTuple, term_descriptor> three_terms;
        std::map<IndexTuple, term_descriptor> two_terms;

    };
}

#endif
