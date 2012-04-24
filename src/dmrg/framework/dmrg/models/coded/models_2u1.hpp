/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_HAMILTONIANS_H
#define APP_HAMILTONIANS_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** FERMI HUBBARD */
template<class Matrix>
class TwoU1_FermiHubbard : public Model<Matrix, TwoU1>
{
public:
    typedef Hamiltonian<Matrix, TwoU1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    
    TwoU1_FermiHubbard(const Lattice& lat, BaseParameters & parms)
    {
        op_t create_up, create_down, destroy_up, destroy_down;
        op_t count_up, count_down, doubly_occ;
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
        
        sign_up.insert_block(Matrix(1, 1, 1), A, A);
        sign_up.insert_block(Matrix(1, 1, -1), B, B);
        sign_up.insert_block(Matrix(1, 1, 1), C, C);
        sign_up.insert_block(Matrix(1, 1, -1), D, D);
        
        sign_down.insert_block(Matrix(1, 1, 1), A, A);
        sign_down.insert_block(Matrix(1, 1, 1), B, B);
        sign_down.insert_block(Matrix(1, 1, -1), C, C);
        sign_down.insert_block(Matrix(1, 1, -1), D, D);
        
        op_t tmp;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                hamterm_t term;
                term.fill_operator = ident;
                term.operators.push_back( std::make_pair(p, parms.get<double>("U")*doubly_occ) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                double ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));
                { // t*cdag_up*c_up
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    gemm(sign_up, create_up, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_up) );
                    terms.push_back(term);
                }
                { // t*c_up*cdag_up
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    gemm(destroy_up, sign_up, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, create_up) );
                    terms.push_back(term);
                }
                { // t*cdag_down*c_down
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    gemm(sign_down, create_down, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_down) );
                    terms.push_back(term);
                }
                { // t*c_down*cdag_down
                    hamterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    gemm(destroy_down, sign_down, tmp); // Note inverse notation because of notation in operator.
                    term.operators.push_back( std::make_pair(p, -ti*tmp) );
                    term.operators.push_back( std::make_pair(*hopto, create_down) );
                    terms.push_back(term);
                }
            }
        }
    }
    
    
    Hamiltonian<Matrix, TwoU1> H () const
    {
        return ham(phys, ident, terms);
    }
    
    Measurements<Matrix, TwoU1> measurements () const
    {
        return Measurements<Matrix, TwoU1>();
    }
    
private:
    Index<TwoU1> phys;
    op_t ident;
    std::vector<hamterm_t> terms;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms.get<double>(key.str()) : parms.get<double>("t");
    }
};

#endif
