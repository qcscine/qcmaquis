/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef QC_TERMMAKER_H
#define QC_TERMMAKER_H

template <class M, class S>
struct TermMaker {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef typename term_descriptor::value_type pos_op_t;
    typedef typename S::subcharge sc_t;

    static bool compare_tag(const pos_op_t& p1, const pos_op_t& p2)
    {
        return p1.first < p2.first;
    }

    static term_descriptor two_term(bool sign, std::vector<tag_type> const & fill_op, value_type scale, pos_t i, pos_t j,
                                     std::vector<tag_type> const & op1, std::vector<tag_type> const & op2,
                                     std::shared_ptr<TagHandler<M, S> > op_table,
                                     Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;
        term.push_back(std::make_pair(i, op1[lat.get_prop<sc_t>("type", i)]));
        term.push_back(std::make_pair(j, op2[lat.get_prop<sc_t>("type", j)]));
        return term;
    }

    static term_descriptor positional_two_term(bool sign, std::vector<tag_type> const & fill_op, value_type scale, pos_t i, pos_t j,
                                     std::vector<tag_type> const & op1, std::vector<tag_type> const & op2,
                                     std::shared_ptr<TagHandler<M, S> > op_table,
                                     Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;

        std::pair<tag_type, value_type> ptag;
        if (i < j) {
            ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", i)], op1[lat.get_prop<sc_t>("type", i)]);
            term.push_back(std::make_pair(i, ptag.first));
            term.push_back(std::make_pair(j, op2[lat.get_prop<sc_t>("type", j)]));
            term.coeff *= ptag.second;
        }
        else {
            ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", j)], op2[lat.get_prop<sc_t>("type", j)]);
            term.push_back(std::make_pair(i, op1[lat.get_prop<sc_t>("type", i)]));
            term.push_back(std::make_pair(j, ptag.first));
            term.coeff *= -ptag.second;
        }
        return term;
    }

    // same, but multiply first two operators
    static term_descriptor positional_two_term(bool sign, std::vector<tag_type> const & fill_op, value_type scale, pos_t i, pos_t j,
                                     std::vector<tag_type> const & op1, std::vector<tag_type> const & op2, std::vector<tag_type> const & op3,
                                     std::shared_ptr<TagHandler<M, S> > op_table,
                                     Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;

        std::pair<tag_type, value_type> pre_ptag;
        pre_ptag = op_table->get_product_tag(op1[lat.get_prop<sc_t>("type",i)], op2[lat.get_prop<sc_t>("type",i)]);

        std::pair<tag_type, value_type> ptag;
        if (i < j) {
            ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", i)], pre_ptag.first);
            term.push_back(std::make_pair(i, ptag.first));
            term.push_back(std::make_pair(j, op3[lat.get_prop<sc_t>("type", j)]));
            term.coeff *= ptag.second * pre_ptag.second;
        }
        else {
            ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", j)], op3[lat.get_prop<sc_t>("type", j)]);
            term.push_back(std::make_pair(i, pre_ptag.first));
            term.push_back(std::make_pair(j, ptag.first));
            term.coeff *= -ptag.second * pre_ptag.second;
        }
        return term;
    }

    static term_descriptor three_term(std::vector<tag_type> const & ident, std::vector<tag_type> const & fill_op,
                                     value_type scale, pos_t pb, pos_t p1, pos_t p2,
                                     std::vector<tag_type> const & opb1, std::vector<tag_type> const & opb2,
                                     std::vector<tag_type> const & ops1, std::vector<tag_type> const & ops2,
                                     std::shared_ptr<TagHandler<M, S> > op_table,
                                     Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = true;
        term.coeff = scale;

        tag_type boson_op;
        tag_type op1 = ops1[lat.get_prop<sc_t>("type", p1)];
        tag_type op2 = ops2[lat.get_prop<sc_t>("type", p2)];
        std::pair<tag_type, value_type> ptag1, ptag2;

        if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
            // if the bosonic operator is in between
            // the fermionic operators, multiply with fill
            ptag1 = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", pb)], opb2[lat.get_prop<sc_t>("type", pb)]);
            term.coeff *= ptag1.second;
            ptag2 = op_table->get_product_tag(ptag1.first, opb1[lat.get_prop<sc_t>("type", pb)]);
            term.coeff *= ptag2.second;
            boson_op = ptag2.first;
        }
        else {
            ptag1 = op_table->get_product_tag(opb2[lat.get_prop<sc_t>("type", pb)], opb1[lat.get_prop<sc_t>("type", pb)]);
            boson_op = ptag1.first;
            term.coeff *= ptag1.second;
        }

        if (p1 < p2) {
            ptag1 = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", p1)], ops1[lat.get_prop<sc_t>("type", p1)]);
            op1 = ptag1.first;
            term.coeff *= ptag1.second;
        }
        else {
            ptag1 = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", p2)], ops2[lat.get_prop<sc_t>("type", p2)]);
            op2 = ptag1.first;
            term.coeff *= -ptag1.second;
        }

        std::vector<pos_op_t> sterm;
        sterm.push_back( std::make_pair(pb, boson_op) );
        sterm.push_back( std::make_pair(p1, op1) );
        sterm.push_back( std::make_pair(p2, op2) );
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        term.push_back(sterm[0]);
        term.push_back(sterm[1]);
        term.push_back(sterm[2]);

        return term;
    }

    static term_descriptor four_term(std::vector<tag_type> const & ident, std::vector<tag_type> const & fill_op,
                                value_type scale, pos_t i, pos_t j, pos_t k, pos_t l,
                                std::vector<tag_type> const & op_i, std::vector<tag_type> const & op_j,
                                std::vector<tag_type> const & op_k, std::vector<tag_type> const & op_l,
                                std::shared_ptr<TagHandler<M, S> > op_table,
                                Lattice const & lat)
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
        sterm.push_back(std::make_pair(i, op_i[lat.get_prop<sc_t>("type", i)]));
        sterm.push_back(std::make_pair(j, op_j[lat.get_prop<sc_t>("type", j)]));
        sterm.push_back(std::make_pair(k, op_k[lat.get_prop<sc_t>("type", k)]));
        sterm.push_back(std::make_pair(l, op_l[lat.get_prop<sc_t>("type", l)]));
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        std::pair<tag_type, value_type> ptag;
        ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", sterm[0].first)], sterm[0].second);
        sterm[0].second = ptag.first;
        term.coeff *= ptag.second;
        ptag = op_table->get_product_tag(fill_op[lat.get_prop<sc_t>("type", sterm[2].first)], sterm[2].second);
        sterm[2].second = ptag.first;
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
