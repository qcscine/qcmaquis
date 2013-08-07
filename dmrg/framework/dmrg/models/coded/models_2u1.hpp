/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MODELS_CODED_2U1_H
#define MODELS_CODED_2U1_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** FERMI HUBBARD */
template<class Matrix>
class FermiHubbardTwoU1 : public Model<Matrix, TwoU1>
{
public:
    typedef Hamiltonian<Matrix, TwoU1> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::hamtagterm_t hamtagterm_t;        
    typedef typename ham::op_t op_t;
    typedef typename OPTable<Matrix, TwoU1>::tag_type tag_type;
    typedef typename Matrix::value_type value_type;
    
    FermiHubbardTwoU1(const Lattice& lat_, BaseParameters & parms_) : lat(lat_), parms(parms_)
    {
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
    }
    
    Index<TwoU1> get_phys() const
    {
        return phys;
    }
    
    Hamiltonian<Matrix, TwoU1> H () const
    {
        std::vector<hamtagterm_t> terms;

        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
             count_up_op, count_down_op, doubly_occ_op,
             sign_up_op, sign_down_op, ident_op;

        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        
        ident_op.insert_block(Matrix(1, 1, 1), A, A);
        ident_op.insert_block(Matrix(1, 1, 1), B, B);
        ident_op.insert_block(Matrix(1, 1, 1), C, C);
        ident_op.insert_block(Matrix(1, 1, 1), D, D);
        
        create_up_op.insert_block(Matrix(1, 1, 1), A, B);
        create_up_op.insert_block(Matrix(1, 1, 1), C, D);
        create_down_op.insert_block(Matrix(1, 1, 1), A, C);
        create_down_op.insert_block(Matrix(1, 1, 1), B, D);
        
        destroy_up_op.insert_block(Matrix(1, 1, 1), B, A);
        destroy_up_op.insert_block(Matrix(1, 1, 1), D, C);
        destroy_down_op.insert_block(Matrix(1, 1, 1), C, A);
        destroy_down_op.insert_block(Matrix(1, 1, 1), D, B);
        
        count_up_op.insert_block(Matrix(1, 1, 1), B, B);
        count_up_op.insert_block(Matrix(1, 1, 1), D, D);
        count_down_op.insert_block(Matrix(1, 1, 1), C, C);
        count_down_op.insert_block(Matrix(1, 1, 1), D, D);
        
        doubly_occ_op.insert_block(Matrix(1, 1, 1), D, D);
        
        sign_up_op.insert_block(Matrix(1, 1, 1), A, A);
        sign_up_op.insert_block(Matrix(1, 1, -1), B, B);
        sign_up_op.insert_block(Matrix(1, 1, 1), C, C);
        sign_up_op.insert_block(Matrix(1, 1, -1), D, D);
        
        sign_down_op.insert_block(Matrix(1, 1, 1), A, A);
        sign_down_op.insert_block(Matrix(1, 1, 1), B, B);
        sign_down_op.insert_block(Matrix(1, 1, -1), C, C);
        sign_down_op.insert_block(Matrix(1, 1, -1), D, D);

        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/

        boost::shared_ptr<OPTable<Matrix, TwoU1> > operator_table(new OPTable<Matrix, TwoU1>());
        tag_type create_up, create_down, destroy_up, destroy_down,
              count_up, count_down, doubly_occ,
              ident, sign_up, sign_down;

        #define REGISTER(op) op = operator_table->register_op(op ## _op);

        REGISTER(ident)
        REGISTER(sign_up)
        REGISTER(sign_down)
        REGISTER(create_up)
        REGISTER(create_down)
        REGISTER(destroy_up)
        REGISTER(destroy_down)
        REGISTER(count_up)
        REGISTER(count_down)
        REGISTER(doubly_occ)

        #undef REGISTER
        /**********************************************************************/
        
        value_type U = parms["U"];
        std::pair<tag_type, value_type> ptag;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                hamtagterm_t term;
                term.fill_operator = ident;
                term.scale = U;
                term.operators.push_back( std::make_pair(p, doubly_occ) );
                terms.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                value_type ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));
                { // t*cdag_up*c_up
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    term.scale = -ti;

                    //gemm(sign_up, create_up, tmp); 
                    ptag = operator_table->get_product_tag(sign_up, create_up); // Note inverse notation because of notation in operator.
                    term.scale *= ptag.second;

                    term.operators.push_back( std::make_pair(p, ptag.first) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_up) );
                    terms.push_back(term);
                }
                { // t*c_up*cdag_up
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_up;
                    term.scale = -ti;

                    //gemm(destroy_up, sign_up, tmp);
                    ptag = operator_table->get_product_tag(destroy_up, sign_up); // Note inverse notation because of notation in operator.
                    term.scale *= ptag.second;

                    term.operators.push_back( std::make_pair(p, ptag.first) );
                    term.operators.push_back( std::make_pair(*hopto, create_up) );
                    terms.push_back(term);
                }
                { // t*cdag_down*c_down
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    term.scale = -ti;

                    //gemm(sign_down, create_down, tmp);
                    ptag = operator_table->get_product_tag(sign_down, create_down); // Note inverse notation because of notation in operator.
                    term.scale *= ptag.second;

                    term.operators.push_back( std::make_pair(p, ptag.first) );
                    term.operators.push_back( std::make_pair(*hopto, destroy_down) );
                    terms.push_back(term);
                }
                { // t*c_down*cdag_down
                    hamtagterm_t term;
                    term.with_sign = true;
                    term.fill_operator = sign_down;
                    term.scale = -ti;

                    //gemm(destroy_down, sign_down, tmp);
                    ptag = operator_table->get_product_tag(destroy_down, sign_down); // Note inverse notation because of notation in operator.
                    term.scale *= ptag.second;

                    term.operators.push_back( std::make_pair(p, ptag.first) );
                    term.operators.push_back( std::make_pair(*hopto, create_down) );
                    terms.push_back(term);
                }
            }
        }

        std::vector<hamterm_t> terms_ops;
        return ham(phys, ident_op, terms_ops, ident, terms, operator_table);
    }
    
    Measurements<Matrix, TwoU1> measurements () const
    {
        return Measurements<Matrix, TwoU1>();
    }
    
private:
    Index<TwoU1> phys;

    Lattice const & lat;
    BaseParameters & parms;
    
    double get_t (BaseParameters & parms, int i) const
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};

#endif
