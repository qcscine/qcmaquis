/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2012 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_HAMILTONIANS_H
#define QC_HAMILTONIANS_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

template<class Matrix>
class qc_model : public Model<Matrix, TwoU1>
{
public:
    typedef Hamiltonian<Matrix, TwoU1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    
    qc_model(const Lattice& lat, BaseParameters & parms)
    {
        op_t sign_up, sign_down;
        
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
        
        ident.insert_block(Matrix(1, 1, 1), A, A);
        ident.insert_block(Matrix(1, 1, 1), B, B);
        ident.insert_block(Matrix(1, 1, 1), C, C);
        ident.insert_block(Matrix(1, 1, 1), D, D);
        
        create_up.insert_block(Matrix(1, 1, 1), A, B);
        create_up.insert_block(Matrix(1, 1, 1), C, D);
        create_down.insert_block(Matrix(1, 1, 1), A, C);
        create_down.insert_block(Matrix(1, 1, 1), B, D);
        
        destroy_up.insert_block(Matrix(1, 1, 1), B, A);
        destroy_up.insert_block(Matrix(1, 1, 1), D, C);
        destroy_down.insert_block(Matrix(1, 1, 1), C, A);
        destroy_down.insert_block(Matrix(1, 1, 1), D, B);
        
        count_up.insert_block(Matrix(1, 1, 1), B, B);
        count_up.insert_block(Matrix(1, 1, 1), D, D);
        count_down.insert_block(Matrix(1, 1, 1), C, C);
        count_down.insert_block(Matrix(1, 1, 1), D, D);
        
        doubly_occ.insert_block(Matrix(1, 1, 1), D, D);
        empty2doubly_occ.insert_block(Matrix(1, 1, 1), A, D);
        doubly_occ2empty.insert_block(Matrix(1, 1, 1), D, A);
        
        sign_up.insert_block(Matrix(1, 1, 1), A, A);
        sign_up.insert_block(Matrix(1, 1, -1), B, B);
        sign_up.insert_block(Matrix(1, 1, 1), C, C);
        sign_up.insert_block(Matrix(1, 1, -1), D, D);
        
        sign_down.insert_block(Matrix(1, 1, 1), A, A);
        sign_down.insert_block(Matrix(1, 1, 1), B, B);
        sign_down.insert_block(Matrix(1, 1, -1), C, C);
        sign_down.insert_block(Matrix(1, 1, -1), D, D);

        gemm(sign_up, sign_down, fill);

        op_t tmp;

        gemm(fill, create_down, tmp);
        create_down = tmp;
        gemm(destroy_down, fill, tmp);
        destroy_down = tmp;

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
                            
    Hamiltonian<Matrix, TwoU1> H () const;
    
    Measurements<Matrix, TwoU1> measurements () const
    {
        return Measurements<Matrix, TwoU1>();
    }
    
private:

    Index<TwoU1> phys;

    op_t create_up, create_down, destroy_up, destroy_down;
    op_t count_up, count_down, doubly_occ, empty2doubly_occ, doubly_occ2empty;
    op_t ident, fill;

    std::vector<double> matrix_elements;
    std::vector<std::vector<int> > idx;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms.get<double>(key.str()) : parms.get<double>("t");
    }

    hamterm_t make_two_term(bool sign, op_t fill_op, int i, int j,
                            op_t op1, op_t op2) const
    {
        hamterm_t term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        term.operators.push_back( std::make_pair(i, op1) );
        term.operators.push_back( std::make_pair(j, op2) );
        return term;
    }

    hamterm_t make_positional_two_term(bool sign, op_t fill_op, int i, int j,
                                       op_t op1, op_t op2) const
    {
        hamterm_t term;
        term.with_sign = sign;
        term.fill_operator = fill_op;

        op_t tmp;
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

    static bool compare(std::pair<size_t, op_t> const & p1, std::pair<size_t, op_t> const & p2) {
        return p1.first < p2.first;
    }

    hamterm_t make_three_term(bool sign, op_t fill_op, int pb, int p1, int p2,
                              op_t opb1, op_t opb2, op_t op1, op_t op2) const
    {
        hamterm_t term;
        term.with_sign = sign;
        term.fill_operator = fill_op;
        op_t tmp, boson_op;

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

        std::vector<std::pair<size_t, op_t> > sterm;
        sterm.push_back( std::make_pair(pb, boson_op) );
        sterm.push_back( std::make_pair(p1, op1) );
        sterm.push_back( std::make_pair(p2, op2) );
        std::sort(sterm.begin(), sterm.end(), compare);

        term.operators.push_back(sterm[0]);
        if (pb == sterm[0].first)
            for(int ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, ident) );
        else
            for(int ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, fill) );

        term.operators.push_back(sterm[1]);
        if (pb == sterm[2].first)
            for(int ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, ident) );
        else 
            for(int ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
                term.operators.push_back( std::make_pair(ipad, fill) );

        term.operators.push_back(sterm[2]);

        return term;
    }

    hamterm_t make_positional_four_term(bool sign, op_t fill_op, int i, int j, int k, int l,
                  op_t op_i, op_t op_j, op_t op_k, op_t op_l) const
    {
        hamterm_t term;
        term.with_sign = sign;
        term.fill_operator = fill_op ;
        op_t tmp;

        // Simple O(n^2) algorithm to determine sign of permutation
        int idx[] = { i,j,k,l };
        int inv_count=0, n=4;
        for(int c1 = 0; c1 < n - 1; c1++)
            for(int c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        std::vector<std::pair<size_t, op_t> > sterm;
        sterm.push_back(std::make_pair(i, op_i));
        sterm.push_back(std::make_pair(j, op_j));
        sterm.push_back(std::make_pair(k, op_k));
        sterm.push_back(std::make_pair(l, op_l));
        std::sort(sterm.begin(), sterm.end(), compare);

        gemm(fill, sterm[0].second, tmp);
        sterm[0].second = tmp;
        gemm(fill, sterm[2].second, tmp);
        sterm[2].second = tmp;
        
        if (inv_count % 2) {
            sterm[0].second *= -1;
            term.operators.push_back(sterm[0]);
        }
        else
            term.operators.push_back(sterm[0]);
        for(int ipad=sterm[0].first +1; ipad < sterm[1].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, fill) );

        term.operators.push_back(sterm[1]);
        for(int ipad=sterm[1].first +1; ipad < sterm[2].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, ident) );

        term.operators.push_back(sterm[2]);
        for(int ipad=sterm[2].first +1; ipad < sterm[3].first; ++ipad)
            term.operators.push_back( std::make_pair(ipad, fill) );

        term.operators.push_back(sterm[3]);
        return term;
    }
};

#include "model_qc.hpp"

#endif
