/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_TERMMAKER_H
#define QC_TERMMAKER_H

template <class M>
struct TermMaker {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;

    struct hamtagterm_t { typedef typename Hamiltonian<M, TwoU1>::hamtagterm_t type; };

    struct op_tag_t { typedef typename OPTagTable<M, TwoU1>::op_tag_t type; };

    struct op_tag_pair_t { typedef typename generate_mpo::Operator_Tag_Term<M, TwoU1>::op_pair_t type; };

    static bool compare_tag(typename op_tag_pair_t::type const & p1,
                            typename op_tag_pair_t::type const & p2)
    {
        return p1.first < p2.first;
    }

    static typename hamtagterm_t::type two_term(bool sign, typename op_tag_t::type fill_op, value_type scale, pos_t i, pos_t j,
                                     typename op_tag_t::type op1, typename op_tag_t::type op2,
                                     boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
    {
        typename hamtagterm_t::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        term.scale = scale;
        term.operators.push_back(std::make_pair(i, op1));
        term.operators.push_back(std::make_pair(j, op2));
        return term;
    }

    static typename hamtagterm_t::type positional_two_term(bool sign, typename op_tag_t::type fill_op, value_type scale, pos_t i, pos_t j,
                                     typename op_tag_t::type op1, typename op_tag_t::type op2,
                                     boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
    {
        typename hamtagterm_t::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        term.scale = scale;

        //typename op_t::type tmp;
        std::pair<typename op_tag_t::type, value_type> ptag;
        if (i < j) {
            ptag = op_table->get_product_tag(fill_op, op1);
            term.operators.push_back(std::make_pair(i, ptag.first));
            term.operators.push_back(std::make_pair(j, op2));
            term.scale *= ptag.second;
        }
        else {
            ptag = op_table->get_product_tag(fill_op, op2);
            term.operators.push_back(std::make_pair(i, op1));
            term.operators.push_back(std::make_pair(j, ptag.first));
            term.scale *= -ptag.second;
        }
        return term;
    }

    static typename hamtagterm_t::type three_term(typename op_tag_t::type ident, typename op_tag_t::type fill_op,
                                     value_type scale, pos_t pb, pos_t p1, pos_t p2,
                                     typename op_tag_t::type opb1, typename op_tag_t::type opb2,
                                     typename op_tag_t::type op1,  typename op_tag_t::type op2,
                                     boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
    {
        typename hamtagterm_t::type term;
        term.with_sign = true;
        term.fill_operator = fill_op;
        term.scale = scale;

        typename op_tag_t::type tmp, boson_op;
        std::pair<typename op_tag_t::type, value_type> ptag1, ptag2;

        if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
            // if the bosonic operator is in between
            // the fermionic operators, multiply with fill
            ptag1 = op_table->get_product_tag(fill_op, opb2);
            term.scale *= ptag1.second;
            ptag2 = op_table->get_product_tag(ptag1.first, opb1);
            term.scale *= ptag2.second;
            boson_op = ptag2.first;
        }
        else {
            ptag1 = op_table->get_product_tag(opb2, opb1);
            boson_op = ptag1.first;
            term.scale *= ptag1.second; 
        }
        
        if (p1 < p2) {
            ptag1 = op_table->get_product_tag(fill_op, op1); 
            op1 = ptag1.first;
            term.scale *= ptag1.second;
        }
        else {
            ptag1 = op_table->get_product_tag(fill_op, op2); 
            op2 = ptag1.first;
            term.scale *= -ptag1.second;
        }

        std::vector<typename op_tag_pair_t::type> sterm;
        sterm.push_back( std::make_pair(pb, boson_op) );
        sterm.push_back( std::make_pair(p1, op1) );
        sterm.push_back( std::make_pair(p2, op2) );
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        term.operators.push_back(sterm[0]);
        if (pb == sterm[0].first)
            for(pos_t ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, ident) );
        else
            for(pos_t ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, fill_op) );

        term.operators.push_back(sterm[1]);
        if (pb == sterm[2].first)
            for(pos_t ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, ident) );
        else 
            for(pos_t ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, fill_op) );

        term.operators.push_back(sterm[2]);

        return term;
    }

    static typename hamtagterm_t::type four_term(typename op_tag_t::type ident, typename op_tag_t::type fill_op,
                                value_type scale, pos_t i, pos_t j, pos_t k, pos_t l,
                                typename op_tag_t::type op_i, typename op_tag_t::type op_j,
                                typename op_tag_t::type op_k, typename op_tag_t::type op_l,
                                boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
    {
        typename hamtagterm_t::type term;
        term.with_sign = true;
        term.fill_operator = fill_op;
        term.scale = scale;

        // Simple O(n^2) algorithm to determine sign of permutation
        pos_t idx[] = { i,j,k,l };
        pos_t inv_count=0, n=4;
        for(pos_t c1 = 0; c1 < n - 1; c1++)
            for(pos_t c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        std::vector<typename op_tag_pair_t::type> sterm;
        sterm.push_back(std::make_pair(i, op_i));
        sterm.push_back(std::make_pair(j, op_j));
        sterm.push_back(std::make_pair(k, op_k));
        sterm.push_back(std::make_pair(l, op_l));
        std::sort(sterm.begin(), sterm.end(), compare_tag);

        std::pair<typename op_tag_t::type, value_type> ptag;
        ptag = op_table->get_product_tag(fill_op, sterm[0].second);
        sterm[0].second = ptag.first;
        term.scale *= ptag.second;
        ptag = op_table->get_product_tag(fill_op, sterm[2].second);
        sterm[2].second = ptag.first;
        term.scale *= ptag.second;
        
        if (inv_count % 2)
            term.scale = -term.scale;

        term.operators.push_back(sterm[0]);

        for(pos_t ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, fill_op) );

        term.operators.push_back(sterm[1]);
        for(pos_t ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, ident) );

        term.operators.push_back(sterm[2]);
        for(pos_t ipad=sterm[2].first +1; ipad < sterm[3].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, fill_op) );

        term.operators.push_back(sterm[3]);
        return term;
    }
};

#endif
