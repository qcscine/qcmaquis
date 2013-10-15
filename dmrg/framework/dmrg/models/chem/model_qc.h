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
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/term_maker.h"
#include "dmrg/models/chem/chem_detail.h"
#include "dmrg/models/chem/irrep_manager.h"


template<class Matrix, class SymmGroup>
class qc_model : public Model<Matrix, SymmGroup>, HamiltonianTraits
{
    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:
    
    qc_model(Lattice const & lat_, BaseParameters & parms_) : lat(lat_), parms(parms_)
    {
        typename SymmGroup::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
        
    }

    Index<SymmGroup> get_phys() const
    {
        return phys;
    }
                            
    Hamiltonian<Matrix, SymmGroup> H () const
    {
        return H_impl<Matrix>();
    }

    /* Disabled - need to implement iterators for one_matrix */
    /*
    Hamiltonian<one_matrix, SymmGroup> H_chem () const
    {
        return H_impl<one_matrix>();
    }
    */
    
    Measurements<Matrix, SymmGroup> measurements () const
    {
        typedef typename Hamiltonian<Matrix, SymmGroup>::op_t op_t;
        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
             count_up_op, count_down_op, docc_op, e2d_op, d2e_op,
             swap_d2u_op, swap_u2d_op,
             create_up_count_down_op, create_down_count_up_op, destroy_up_count_down_op, destroy_down_count_up_op,
             ident_op, fill_op;

        typename SymmGroup::charge A(0), B(0), C(0), D(1);
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

        docc_op.insert_block(Matrix(1, 1, 1), D, D);

        e2d_op.insert_block(Matrix(1, 1, 1), A, D);
        d2e_op.insert_block(Matrix(1, 1, 1), D, A);

        fill_op.insert_block(Matrix(1, 1, 1), A, A);
        fill_op.insert_block(Matrix(1, 1, -1), B, B);
        fill_op.insert_block(Matrix(1, 1, -1), C, C);
        fill_op.insert_block(Matrix(1, 1, 1), D, D);

        op_t tmp;

        gemm(fill_op, create_down_op, tmp);
        create_down_op = tmp;
        gemm(destroy_down_op, fill_op, tmp);
        destroy_down_op = tmp;

        gemm(create_up_op, destroy_down_op, swap_d2u_op);
        gemm(destroy_up_op, create_down_op, swap_u2d_op);
        gemm(count_down_op, create_up_op, create_up_count_down_op);
        gemm(count_up_op, create_down_op, create_down_count_up_op);
        gemm(count_down_op, destroy_up_op, destroy_up_count_down_op);
        gemm(count_up_op, destroy_down_op, destroy_down_count_up_op);

        typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
        Measurements<Matrix, SymmGroup> meas;
        meas.set_identity(ident_op);

        {
            boost::regex expression("^MEASURE_LOCAL\\[(.*)]$");
            boost::smatch what;
            for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
                std::string lhs = it->key();
                if (boost::regex_match(lhs, what, expression)) {

                    mterm_t term;
                    term.type = mterm_t::Local;
                    term.name = what.str(1);

                    if (it->value() == "Nup")
                        term.operators.push_back(std::make_pair(count_up_op, false));
                    else if (it->value() == "Ndown")
                        term.operators.push_back(std::make_pair(count_down_op, false));
                    else if (it->value() == "Nup*Ndown" || it->value() == "docc")
                        term.operators.push_back(std::make_pair(docc_op, false));
                    else
                        throw std::runtime_error("Invalid observable\nLocal measurements supported so far are \"Nup\" and \"Ndown\"\n");

                    meas.add_term(term);
                }
            }
        }

        {
        boost::regex expression("^MEASURE_CORRELATIONS\\[(.*)]$");
        boost::regex expression_half("^MEASURE_HALF_CORRELATIONS\\[(.*)]$");
        boost::regex expression_nn("^MEASURE_NN_CORRELATIONS\\[(.*)]$");
        boost::regex expression_halfnn("^MEASURE_HALF_NN_CORRELATIONS\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            std::string value;

            mterm_t term;
            term.fill_operator = ident_op;

            bool found = false;
            if (boost::regex_match(lhs, what, expression)) {
                value = it->value();
                found = true;
                term.name = what.str(1);
                term.type = mterm_t::Correlation;
            }
            if (boost::regex_match(lhs, what, expression_half)) {
                value = it->value();
                found = true;
                term.name = what.str(1);
                term.type = mterm_t::HalfCorrelation;
            }
            if (boost::regex_match(lhs, what, expression_nn)) {
                value = it->value();
                found = true;
                term.name = what.str(1);
                term.type = mterm_t::CorrelationNN;
            }
            if (boost::regex_match(lhs, what, expression_halfnn)) {
                value = it->value();
                found = true;
                term.name = what.str(1);
                term.type = mterm_t::HalfCorrelationNN;
            }
            if (found) {

                int f_ops = 0;

                /// split op1:op2:...@p1,p2,p3,... into {op1:op2:...}, {p1,p2,p3,...}
                std::vector<std::string> value_split;
                boost::split( value_split, value, boost::is_any_of("@"));

                /// parse operators op1:op2:...
                boost::char_separator<char> sep(":");
                tokenizer corr_tokens(value_split[0], sep);
                for (tokenizer::iterator it2=corr_tokens.begin();
                     it2 != corr_tokens.end();
                     it2++)
                {
                    if (*it2 == "c_up") {
                        term.operators.push_back( std::make_pair(destroy_up_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_down") {
                        term.operators.push_back( std::make_pair(destroy_down_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_up") {
                        term.operators.push_back( std::make_pair(create_up_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_down") {
                        term.operators.push_back( std::make_pair(create_down_op, true) );
                        ++f_ops;
                    }

                    else if (*it2 == "Nup") {
                        term.operators.push_back( std::make_pair(count_up_op, false) );
                    }
                    else if (*it2 == "Ndown") {
                        term.operators.push_back( std::make_pair(count_down_op, false) );
                    }
                    else if (*it2 == "docc" || *it2 == "Nup*Ndown") {
                        term.operators.push_back( std::make_pair(docc_op, false) );
                    }
                    else if (*it2 == "cdag_up*c_down") {
                        term.operators.push_back( std::make_pair(swap_d2u_op, false) );
                    }
                    else if (*it2 == "cdag_down*c_up") {
                        term.operators.push_back( std::make_pair(swap_u2d_op, false) );
                    }

                    else if (*it2 == "cdag_up*cdag_down") {
                        term.operators.push_back( std::make_pair(e2d_op, false) );
                    }
                    else if (*it2 == "c_up*c_down") {
                        term.operators.push_back( std::make_pair(d2e_op, false) );
                    }

                    else if (*it2 == "cdag_up*Ndown") {
                        term.operators.push_back( std::make_pair(create_up_count_down_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_down*Nup") {
                        term.operators.push_back( std::make_pair(create_down_count_up_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_up*Ndown") {
                        term.operators.push_back( std::make_pair(destroy_up_count_down_op, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_down*Nup") {
                        term.operators.push_back( std::make_pair(destroy_down_count_up_op, true) );
                        ++f_ops;
                    }
                    else
                        throw std::runtime_error("Unrecognized operator in correlation measurement: " 
                                                    + boost::lexical_cast<std::string>(*it2) + "\n");

                }
                //if (term.operators.size() == 1) {
                //    term.operators.push_back(term.operators[0]);
                //    if (term.operators[1].second) ++f_ops;
                //}


                if (f_ops > 0)
                    term.fill_operator = fill_op;

                if (f_ops % 2 != 0)
                    throw std::runtime_error("In " + term.name + ": Number of fermionic operators has to be even in correlation measurements.");

                /// parse positions p1,p2,p3,... (or `space`)
                if (value_split.size() > 1) {
                    boost::char_separator<char> pos_sep(", ");
                    tokenizer pos_tokens(value_split[1], pos_sep);
                    term.positions.resize(1);
                    std::transform(pos_tokens.begin(), pos_tokens.end(), std::back_inserter(term.positions[0]),
                                   static_cast<std::size_t (*)(std::string const&)>(boost::lexical_cast<std::size_t, std::string>));
                }

                meas.add_term(term);
            }
        }
        }
        return meas;
    }

    template <class M>
    Hamiltonian<M, SymmGroup> H_impl () const;

private:
    Index<SymmGroup> phys;

    Lattice const & lat;
    BaseParameters & parms;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }

};


#include "dmrg/models/chem/model_qc.hpp"

#endif
