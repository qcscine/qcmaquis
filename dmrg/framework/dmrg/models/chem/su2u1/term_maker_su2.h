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

#ifndef QC_TERMMAKER_SU2_H
#define QC_TERMMAKER_SU2_H

template <class M, class S>
struct TermMakerSU2 {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef typename term_descriptor::value_type pos_op_t;

    static bool compare_tag(pos_op_t p1,
                            pos_op_t p2)
    {
        return boost::tuples::get<0>(p1) < boost::tuples::get<0>(p2);
    }

    static term_descriptor two_term(bool sign, tag_type full_ident, value_type scale, pos_t i, pos_t j,
                                     tag_type op1, tag_type op2)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.full_identity = full_ident;
        term.coeff = scale;
        term.push_back(boost::make_tuple(i, op1));
        term.push_back(boost::make_tuple(j, op2));
        return term;
    }

    static term_descriptor positional_two_term(bool sign, pos_t lat_size, tag_type ident, tag_type fill, value_type scale, pos_t i, pos_t j,
                                     tag_type op1, tag_type op1_fill, tag_type op2, tag_type op2_fill)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.full_identity = fill;
        term.coeff = scale;

        tag_type op1_use = (i<j) ? op1_fill : op2_fill;
        tag_type op2_use = (i<j) ? op2 : op1;
        if (j<i && sign) term.coeff = -term.coeff;

        int start = std::min(i,j), end = std::max(i,j);
        //for (int fs=0; fs < start; ++fs)
        //    term.push_back( boost::make_tuple(fs, ident) );
        term.push_back( boost::make_tuple(start, op1_use) );

        //for (int fs = start+1; fs < end; ++fs)
        //    term.push_back( boost::make_tuple(fs, fill) );
        term.push_back( boost::make_tuple(end, op2_use) );

        //for (int fs = end+1; fs < lat_size; ++fs)
        //    term.push_back( boost::make_tuple(fs, ident) );

        return term;
    }

    static term_descriptor three_term(tag_type ident, tag_type fill_op,
                                      value_type scale, pos_t pb, pos_t p1, pos_t p2,
                                      tag_type boson_op, tag_type boson_op_fill,
                                      tag_type op1, tag_type op1_fill, tag_type op2, tag_type op2_fill)
    {
        term_descriptor term;
        term.is_fermionic = true;
        term.coeff = scale;
        term.full_identity = ident;

        tag_type boson_op_use, op1_use, op2_use;

        op1_use = (p1<p2) ? op1_fill : op2_fill;
        op2_use = (p1<p2) ? op2 : op1;

        if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
            // if the bosonic operator is in between
            // the fermionic operators, use fill-multiplied version
            boson_op_use = boson_op_fill;
        }
        else {
            boson_op_use = boson_op;
        }

        if (p2<p1) term.coeff = -term.coeff;

        std::vector<pos_op_t> sterm;
        sterm.push_back( boost::make_tuple(pb, boson_op_use) );
        sterm.push_back( boost::make_tuple(std::min(p1,p2), op1_use) );
        sterm.push_back( boost::make_tuple(std::max(p1,p2), op2_use) );
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        term.push_back(sterm[0]);
        term.push_back(sterm[1]);
        term.push_back(sterm[2]);

        return term;
    }

    static term_descriptor four_term(tag_type ident, tag_type fill_op,
                                value_type scale, pos_t i, pos_t j, pos_t k, pos_t l,
                                tag_type op_i, tag_type op_j,
                                tag_type op_k, tag_type op_l,
                                boost::shared_ptr<TagHandler<M, S> > op_table)
    {
        term_descriptor term;
        term.is_fermionic = true;
        term.coeff = scale;

        // Simple O(n^2) algorithm to determine sign of permutation
        pos_t idx[] = { i,j,k,l };
        pos_t inv_count=0, n=4;
        for(pos_t c1 = 0; c1 < n - 1; c1++)
            for(pos_t c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        std::vector<pos_op_t> sterm;
        sterm.push_back(boost::make_tuple(i, op_i));
        sterm.push_back(boost::make_tuple(j, op_j));
        sterm.push_back(boost::make_tuple(k, op_k));
        sterm.push_back(boost::make_tuple(l, op_l));
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        std::pair<tag_type, value_type> ptag;
        ptag = op_table->get_product_tag(fill_op, boost::tuples::get<1>(sterm[0]));
        boost::tuples::get<1>(sterm[0]) = ptag.first;
        term.coeff *= ptag.second;
        ptag = op_table->get_product_tag(fill_op, boost::tuples::get<1>(sterm[2]));
        boost::tuples::get<1>(sterm[2]) = ptag.first;
        term.coeff *= ptag.second;
        
        if (inv_count % 2)
            term.coeff = -term.coeff;

        term.push_back(sterm[0]);
        term.push_back(sterm[1]);
        term.push_back(sterm[2]);
        term.push_back(sterm[3]);
        return term;
    }
};

#endif
