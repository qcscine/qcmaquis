/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_HAMILTONIANS_H
#define QC_HAMILTONIANS_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <boost/shared_ptr.hpp>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"


template<class Matrix>
class qc_model : public Model<Matrix, TwoU1>
{
    typedef typename Lattice::pos_t pos_t;
    
    template <class M>
    struct hamterm_t { typedef typename Hamiltonian<M, TwoU1>::hamterm_t type; };

    template <class M>
    struct hamtagterm_t { typedef typename Hamiltonian<M, TwoU1>::hamtagterm_t type; };

    template <class M>
    struct op_t { typedef typename Hamiltonian<M, TwoU1>::op_t type; };

    template <class M>
    struct op_tag_t { typedef typename OPTagTable<M, TwoU1>::op_tag_t type; };

    template <class M>
    struct op_pair_t { typedef typename generate_mpo::Operator_Term<M, TwoU1>::op_pair_t type; };

    template <class M>
    struct op_tag_pair_t { typedef typename generate_mpo::Operator_Tag_Term<M, TwoU1>::op_pair_t type; };

    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:
    
    qc_model(const Lattice& lat, BaseParameters & parms)
    {
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
        
        // ********************************************************************
        // *** Parse orbital data *********************************************
        // ********************************************************************

        std::ifstream orb_file;
        orb_file.open(parms.get<std::string>("integral_file").c_str());
        for(int i=0; i < 4; ++i){
            orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
        }

        std::vector<double> raw;
        std::copy(std::istream_iterator<double>(orb_file), std::istream_iterator<double>(),
                    std::back_inserter(raw));

        std::vector<double>::iterator it = raw.begin(); 
        while (it != raw.end()) {
            matrix_elements.push_back(*it++);
            idx.push_back(std::vector<int>(it, it + 4));
            it += 4;
        }

#ifndef NDEBUG
        for(std::size_t m=0; m < matrix_elements.size(); ++m)
        {
            assert( *max_element(idx[m].begin(), idx[m].end()) <= lat.size() );
        }   
#endif

    } // ctor

    Index<TwoU1> get_phys() const
    {
        return phys;
    }
                            
    Hamiltonian<Matrix, TwoU1> H () const
    {
        return H_impl<Matrix>();
    }

    /* Disabled - need to implement iterators for one_matrix */
    /*
    Hamiltonian<one_matrix, TwoU1> H_chem () const
    {
        return H_impl<one_matrix>();
    }
    */
    
    Measurements<Matrix, TwoU1> measurements () const
    {
        return Measurements<Matrix, TwoU1>();
    }

    template <class M>
    Hamiltonian<M, TwoU1> H_impl () const;

private:

    Index<TwoU1> phys;

    std::vector<double> matrix_elements;
    std::vector<std::vector<int> > idx;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms.get<double>(key.str()) : parms.get<double>("t");
    }

    template<class M>
    static bool compare(typename op_pair_t<M>::type const & p1,
                            typename op_pair_t<M>::type const & p2) {
        return p1.first < p2.first;
    }

    template<class M>
    static bool compare_tag(typename op_tag_pair_t<M>::type const & p1,
                            typename op_tag_pair_t<M>::type const & p2) {
        return p1.first < p2.first;
    }

    template <class M>
    struct TermMaker {

        static typename hamtagterm_t<M>::type two_term(bool sign, typename op_tag_t<M>::type fill_op, double scale, pos_t i, pos_t j,
                                         typename op_tag_t<M>::type op1, typename op_tag_t<M>::type op2,
                                         boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
        {
            typename hamtagterm_t<M>::type term;
            term.with_sign = sign;
            term.fill_operator = fill_op;
            term.scale = scale;
            term.operators.push_back(std::make_pair(i, op1));
            term.operators.push_back(std::make_pair(j, op2));
            return term;
        }

        static typename hamtagterm_t<M>::type positional_two_term(bool sign, typename op_tag_t<M>::type fill_op, double scale, pos_t i, pos_t j,
                                         typename op_tag_t<M>::type op1, typename op_tag_t<M>::type op2,
                                         boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
        {
            typename hamtagterm_t<M>::type term;
            term.with_sign = sign;
            term.fill_operator = fill_op;
            term.scale = scale;

            typename op_t<M>::type tmp;
            if (i < j) {
                typename op_tag_t<M>::type tmp = op_table->get_prod_tag(fill_op, op1);
                term.operators.push_back(std::make_pair(i, tmp));
                term.operators.push_back(std::make_pair(j, op2));
            }
            else {
                typename op_tag_t<M>::type tmp = op_table->get_prod_tag(fill_op, op2);
                term.operators.push_back(std::make_pair(i, op1));
                term.operators.push_back(std::make_pair(j, tmp));
                term.scale = -term.scale;
            }
            return term;
        }

        static typename hamtagterm_t<M>::type three_term(bool sign, typename op_tag_t<M>::type ident, typename op_tag_t<M>::type fill_op,
                                         double scale, pos_t pb, pos_t p1, pos_t p2,
                                         typename op_tag_t<M>::type opb1, typename op_tag_t<M>::type opb2,
                                         typename op_tag_t<M>::type op1,  typename op_tag_t<M>::type op2,
                                         boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
        {
            typename hamtagterm_t<M>::type term;
            term.with_sign = sign;
            term.fill_operator = fill_op;
            term.scale = scale;

            typename op_tag_t<M>::type tmp, boson_op;

            if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
                // if the bosonic operator is in between
                // the fermionic operators, multiply with fill
                tmp = op_table->get_prod_tag(fill_op, opb2);
                boson_op = op_table->get_prod_tag(tmp, opb1);
            }
            else {
                boson_op = op_table->get_prod_tag(opb2, opb1);
            }
            
            if (p1 < p2) {
                op1 = op_table->get_prod_tag(fill_op, op1);
            }
            else {
                op2 = op_table->get_prod_tag(fill_op, op2);
                term.scale = -term.scale;
            }

            std::vector<typename op_tag_pair_t<M>::type> sterm;
            sterm.push_back( std::make_pair(pb, boson_op) );
            sterm.push_back( std::make_pair(p1, op1) );
            sterm.push_back( std::make_pair(p2, op2) );
            std::sort(sterm.begin(), sterm.end(), compare_tag<M>);

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

        static typename hamtagterm_t<M>::type four_term(bool sign, typename op_tag_t<M>::type ident, typename op_tag_t<M>::type fill_op,
                                    double scale, pos_t i, pos_t j, pos_t k, pos_t l,
                                    typename op_tag_t<M>::type op_i, typename op_tag_t<M>::type op_j,
                                    typename op_tag_t<M>::type op_k, typename op_tag_t<M>::type op_l,
                                    boost::shared_ptr<OPTagTable<M, TwoU1> > op_table)
        {
            typename hamtagterm_t<M>::type term;
            term.with_sign = sign;
            term.fill_operator = fill_op;
            term.scale = scale;

            // Simple O(n^2) algorithm to determine sign of permutation
            pos_t idx[] = { i,j,k,l };
            pos_t inv_count=0, n=4;
            for(pos_t c1 = 0; c1 < n - 1; c1++)
                for(pos_t c2 = c1+1; c2 < n; c2++)
                    if(idx[c1] > idx[c2]) inv_count++;

            std::vector<typename op_tag_pair_t<M>::type> sterm;
            sterm.push_back(std::make_pair(i, op_i));
            sterm.push_back(std::make_pair(j, op_j));
            sterm.push_back(std::make_pair(k, op_k));
            sterm.push_back(std::make_pair(l, op_l));
            std::sort(sterm.begin(), sterm.end(), compare_tag<M>);

            sterm[0].second = op_table->get_prod_tag(fill_op, sterm[0].second);
            sterm[2].second = op_table->get_prod_tag(fill_op, sterm[2].second);
            
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

    /*
    template <class M>
    typename hamterm_t<M>::type make_two_term(bool sign, typename op_t<M>::type fill_op, pos_t i, pos_t j,
                                     typename op_t<M>::type op1, typename op_t<M>::type op2) const
    {
        typename hamterm_t<M>::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        term.operators.push_back( std::make_pair(i, op1) );
        term.operators.push_back( std::make_pair(j, op2) );
        return term;
    }

    template <class M>
    typename hamterm_t<M>::type make_positional_two_term(bool sign, typename op_t<M>::type fill_op, pos_t i, pos_t j,
                                     typename op_t<M>::type op1, typename op_t<M>::type op2) const
    {
        typename hamterm_t<M>::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op;

        typename op_t<M>::type tmp;
        if (i < j) {
            gemm(fill_op, op1, tmp);
            term.operators.push_back( std::make_pair(i, tmp) );
            term.operators.push_back( std::make_pair(j, op2) );
        }
        else {
            gemm(fill_op, op2, tmp);
            term.operators.push_back( std::make_pair(i, op1) );
            term.operators.push_back( std::make_pair(j, -1*tmp) );
        }
        return term;
    }

    template<class M>
    typename hamterm_t<M>::type make_three_term(bool sign, typename op_t<M>::type ident, typename op_t<M>::type fill_op,
                                     pos_t pb, pos_t p1, pos_t p2,
                                     typename op_t<M>::type opb1, typename op_t<M>::type opb2,
                                     typename op_t<M>::type op1,  typename op_t<M>::type op2) const
    {  
        typename hamterm_t<M>::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        typename op_t<M>::type tmp, boson_op;

        if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
            // if the bosonic operator is in between
            // the fermionic operators, multiply with fill
            gemm(fill_op, opb2, tmp);
            gemm(tmp, opb1, boson_op);
        }
        else {
            gemm(opb2, opb1, boson_op);
        }

        if (p1 < p2) {
            gemm(fill_op, op1, tmp);
            op1 = tmp;
        }
        else {
            gemm(fill_op, op2, tmp);
            op2 = -1*tmp;
        }

        std::vector<typename op_pair_t<M>::type> sterm;
        sterm.push_back( std::make_pair(pb, boson_op) );
        sterm.push_back( std::make_pair(p1, op1) );
        sterm.push_back( std::make_pair(p2, op2) );
        std::sort(sterm.begin(), sterm.end(), compare<M>);

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

    template <class M>
    typename hamterm_t<M>::type make_positional_four_term(bool sign, typename op_t<M>::type ident, typename op_t<M>::type fill_op,
                                pos_t i, pos_t j, pos_t k, pos_t l,
                                typename op_t<M>::type op_i, typename op_t<M>::type op_j,
                                typename op_t<M>::type op_k, typename op_t<M>::type op_l) const
    {
        typename hamterm_t<M>::type term;
        term.with_sign = sign;
        term.fill_operator = fill_op ;
        typename op_t<M>::type tmp;

        // Simple O(n^2) algorithm to determine sign of permutation
        pos_t idx[] = { i,j,k,l };
        pos_t inv_count=0, n=4;
        for(pos_t c1 = 0; c1 < n - 1; c1++)
            for(pos_t c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        std::vector<typename op_pair_t<M>::type> sterm;
        sterm.push_back(std::make_pair(i, op_i));
        sterm.push_back(std::make_pair(j, op_j));
        sterm.push_back(std::make_pair(k, op_k));
        sterm.push_back(std::make_pair(l, op_l));
        std::sort(sterm.begin(), sterm.end(), compare<M>);

        gemm(fill_op, sterm[0].second, tmp);
        sterm[0].second = tmp;
        gemm(fill_op, sterm[2].second, tmp);
        sterm[2].second = tmp;

        if (inv_count % 2)
            sterm[0].second *= -1;

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
    */

};


#include "model_qc.hpp"

#endif
