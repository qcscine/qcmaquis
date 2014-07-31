/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/term_maker.h"
#include "dmrg/models/chem/chem_detail.h"
#include "dmrg/models/chem/pg_util.h"

template<class Matrix, class SymmGroup>
class qc_model : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:
    
    qc_model(Lattice const & lat_, BaseParameters & parms_);
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill;
    }

    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return chem_detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "create_up")
            return create_up;
        else if (name == "create_down")
            return create_down;
        else if (name == "destroy_up")
            return destroy_up;
        else if (name == "destroy_down")
            return destroy_down;
        else if (name == "count_up")
            return count_up;
        else if (name == "count_down")
            return count_down;
        else if (name == "e2d")
            return e2d;
        else if (name == "d2e")
            return d2e;
        else if (name == "docc")
            return docc;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    measurements_type measurements () const
    {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
             count_up_op, count_down_op, docc_op, e2d_op, d2e_op,
             swap_d2u_op, swap_u2d_op,
             create_up_count_down_op, create_down_count_up_op, destroy_up_count_down_op, destroy_down_count_up_op,
             ident_op, fill_op;

        ident_op = tag_handler->get_op(ident);
        fill_op = tag_handler->get_op(fill);
        create_up_op = tag_handler->get_op(create_up);
        create_down_op = tag_handler->get_op(create_down);
        destroy_up_op = tag_handler->get_op(destroy_up);
        destroy_down_op = tag_handler->get_op(destroy_down);
        count_up_op = tag_handler->get_op(count_up);
        count_down_op = tag_handler->get_op(count_down);
        e2d_op = tag_handler->get_op(e2d);
        d2e_op = tag_handler->get_op(d2e);
        docc_op = tag_handler->get_op(docc);

        gemm(destroy_down_op, create_up_op, swap_d2u_op); // S_plus
        gemm(destroy_up_op, create_down_op, swap_u2d_op); // S_minus

        gemm(count_down_op, create_up_op, create_up_count_down_op);
        gemm(count_up_op, create_down_op, create_down_count_up_op);
        gemm(count_down_op, destroy_up_op, destroy_up_count_down_op);
        gemm(count_up_op, destroy_down_op, destroy_down_count_up_op);

        #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

        GENERATE_SITE_SPECIFIC(ident_op)
        GENERATE_SITE_SPECIFIC(fill_op)
        GENERATE_SITE_SPECIFIC(create_up_op)
        GENERATE_SITE_SPECIFIC(create_down_op)
        GENERATE_SITE_SPECIFIC(destroy_up_op)
        GENERATE_SITE_SPECIFIC(destroy_down_op)
        GENERATE_SITE_SPECIFIC(count_up_op)
        GENERATE_SITE_SPECIFIC(count_down_op)

        GENERATE_SITE_SPECIFIC(e2d_op)
        GENERATE_SITE_SPECIFIC(d2e_op)
        GENERATE_SITE_SPECIFIC(docc_op)

        GENERATE_SITE_SPECIFIC(swap_d2u_op)
        GENERATE_SITE_SPECIFIC(swap_u2d_op)
        GENERATE_SITE_SPECIFIC(create_up_count_down_op)
        GENERATE_SITE_SPECIFIC(create_down_count_up_op)
        GENERATE_SITE_SPECIFIC(destroy_up_count_down_op)
        GENERATE_SITE_SPECIFIC(destroy_down_count_up_op)

        #undef GENERATE_SITE_SPECIFIC

        
        measurements_type meas;

        typedef std::vector<block_matrix<Matrix, SymmGroup> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        {
            boost::regex expression("^MEASURE_LOCAL\\[(.*)]$");
            boost::smatch what;
            for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
                std::string lhs = it->key();
                if (boost::regex_match(lhs, what, expression)) {

                    op_vec meas_op;
                    if (it->value() == "Nup")
                        meas_op = count_up_ops;
                    else if (it->value() == "Ndown")
                        meas_op = count_down_ops;
                    else if (it->value() == "Nup*Ndown" || it->value() == "docc")
                        meas_op = docc_ops;
                    else
                        throw std::runtime_error("Invalid observable\nLocal measurements supported so far are \"Nup\" and \"Ndown\"\n");

                    meas.push_back( new measurements::local<Matrix, SymmGroup>(what.str(1), lat, ident_ops, fill_ops, meas_op) );
                }
            }
        }

        {
        boost::regex expression("^MEASURE_CORRELATIONS\\[(.*)]$");
        boost::regex expression_half("^MEASURE_HALF_CORRELATIONS\\[(.*)]$");
        boost::regex expression_nn("^MEASURE_NN_CORRELATIONS\\[(.*)]$");
        boost::regex expression_halfnn("^MEASURE_HALF_NN_CORRELATIONS\\[(.*)]$");
        boost::regex expression_twoptdm("^MEASURE_TWOPTDM(.*)$");
        boost::regex expression_transition_twoptdm("^MEASURE_TRANSITION_TWOPTDM(.*)$");
        boost::regex expression_threeptdm("^THREEPTDM(.*)$");
        boost::regex expression_transition_threeptdm("^MEASURE_TRANSITION_THREEPTDM(.*)$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();

            std::string name, value;
            bool half_only, nearest_neighbors_only;
            if (boost::regex_match(lhs, what, expression)) {
                value = it->value();
                name = what.str(1);
                half_only = false;
                nearest_neighbors_only = false;
            }
            if (boost::regex_match(lhs, what, expression_half)) {
                value = it->value();
                name = what.str(1);
                half_only = true;
                nearest_neighbors_only = false;
            }
            if (boost::regex_match(lhs, what, expression_nn)) {
                value = it->value();
                name = what.str(1);
                half_only = false;
                nearest_neighbors_only = true;
            }
            if (boost::regex_match(lhs, what, expression_halfnn)) {
                value = it->value();
                name = what.str(1);
                half_only = true;
                nearest_neighbors_only = true;
            }
            if (boost::regex_match(lhs, what, expression_twoptdm) ||
                    boost::regex_match(lhs, what, expression_transition_twoptdm)) {

                maquis::cout << "matched twoptdm measurement with lhs = " << lhs << std::endl;
                std::string bra_ckp("");
                if(lhs == "MEASURE_TRANSITION_TWOPTDM"){
                    name = "transition_twoptdm";
                    bra_ckp = it->value();
                }
                else
                    name = "twoptdm";
                maquis::cout << "Added 2pdm for measuring: name == "<< name <<"\n";

                std::vector<bond_element> synchronous_meas_operators;
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                half_only = true;
                nearest_neighbors_only = false;
                std::vector<pos_t> positions;
                meas.push_back( new measurements::NRankRDM<Matrix, SymmGroup>(name, lat, ident_ops, fill_ops, synchronous_meas_operators,
                                                                              half_only, nearest_neighbors_only, positions, bra_ckp));
            }
            else if (boost::regex_match(lhs, what, expression_threeptdm) ||
                    boost::regex_match(lhs, what, expression_transition_threeptdm)) {
  
                std::string bra_ckp("");
                if(lhs == "MEASURE_TRANSITION_THREEPTDM"){
                    name = "transition_threeptdm";
                    bra_ckp = it->value();
                }
                else
                    name = "threeptdm";
                maquis::cout << "Added 3pdm for measuring: name == "<< name <<"\n";
  
                std::vector<bond_element> synchronous_meas_operators;
                {
                //* c+a c+a c+a ca ca ca *//
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                //* c+a c+a c+b ca ca cb *//
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                //* c+b c+b c+a cb cb ca *//
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_up_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                //* c+b c+b c+b cb cb cb *//
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(create_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                half_only = true;
                nearest_neighbors_only = false;
                std::vector<pos_t> positions;
                meas.push_back( new measurements::NRankRDM<Matrix, SymmGroup>(name, lat, ident_ops, fill_ops, synchronous_meas_operators,
                                                                              half_only, nearest_neighbors_only, positions, bra_ckp));
            }
            else if (!name.empty()) {

                int f_ops = 0;
                bond_element meas_operators;
                
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
                        meas_operators.push_back( std::make_pair(destroy_up_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_down") {
                        meas_operators.push_back( std::make_pair(destroy_down_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_up") {
                        meas_operators.push_back( std::make_pair(create_up_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_down") {
                        meas_operators.push_back( std::make_pair(create_down_ops, true) );
                        ++f_ops;
                    }

                    else if (*it2 == "id" || *it2 == "Id") {
                        meas_operators.push_back( std::make_pair(ident_ops, false) );
                    }
                    else if (*it2 == "Nup") {
                        meas_operators.push_back( std::make_pair(count_up_ops, false) );
                    }
                    else if (*it2 == "Ndown") {
                        meas_operators.push_back( std::make_pair(count_down_ops, false) );
                    }
                    else if (*it2 == "docc" || *it2 == "Nup*Ndown") {
                        meas_operators.push_back( std::make_pair(docc_ops, false) );
                    }
                    else if (*it2 == "cdag_up*c_down" || *it2 == "splus") {
                        meas_operators.push_back( std::make_pair(swap_d2u_ops, false) );
                    }
                    else if (*it2 == "cdag_down*c_up" || *it2 == "sminus") {
                        meas_operators.push_back( std::make_pair(swap_u2d_ops, false) );
                    }

                    else if (*it2 == "cdag_up*cdag_down") {
                        meas_operators.push_back( std::make_pair(e2d_ops, false) );
                    }
                    else if (*it2 == "c_up*c_down") {
                        meas_operators.push_back( std::make_pair(d2e_ops, false) );
                    }

                    else if (*it2 == "cdag_up*Ndown") {
                        meas_operators.push_back( std::make_pair(create_up_count_down_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_down*Nup") {
                        meas_operators.push_back( std::make_pair(create_down_count_up_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_up*Ndown") {
                        meas_operators.push_back( std::make_pair(destroy_up_count_down_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_down*Nup") {
                        meas_operators.push_back( std::make_pair(destroy_down_count_up_ops, true) );
                        ++f_ops;
                    }
                    else
                        throw std::runtime_error("Unrecognized operator in correlation measurement: " 
                                                    + boost::lexical_cast<std::string>(*it2) + "\n");
                }

                if (f_ops % 2 != 0)
                    throw std::runtime_error("In " + name + ": Number of fermionic operators has to be even in correlation measurements.");

                /// parse positions p1,p2,p3,... (or `space`)
                std::vector<pos_t> positions;
                if (value_split.size() > 1) {
                    boost::char_separator<char> pos_sep(", ");
                    tokenizer pos_tokens(value_split[1], pos_sep);
                    std::transform(pos_tokens.begin(), pos_tokens.end(), std::back_inserter(positions),
                                   static_cast<pos_t (*)(std::string const&)>(boost::lexical_cast<pos_t, std::string>));
                }
                
                std::vector<bond_element> synchronous_meas_operators;
                synchronous_meas_operators.push_back(meas_operators);
                maquis::cout << "adding measurement " << name << std::endl;
                meas.push_back( new measurements::NRankRDM<Matrix, SymmGroup>(name, lat, ident_ops, fill_ops, synchronous_meas_operators,
                                                                              half_only, nearest_neighbors_only, positions));
            }
        }
        }
        return meas;
    }

private:
    Lattice const & lat;
    BaseParameters & parms;
    Index<SymmGroup> phys;

    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    tag_type ident, fill,
             create_up, create_down, destroy_up, destroy_down,
             count_up, count_down, docc, e2d, d2e;

    typename SymmGroup::subcharge max_irrep;

    std::vector<op_t> generate_site_specific_ops(op_t const & op) const
    {
        PGDecorator<SymmGroup> set_symm;
        std::vector<op_t> ret;
        for (typename SymmGroup::subcharge sc=0; sc < max_irrep+1; ++sc) {
            op_t mod(set_symm(op.left_basis(), sc), set_symm(op.right_basis(), sc));
            for (std::size_t b = 0; b < op.n_blocks(); ++b)
                mod[b] = op[b];

            ret.push_back(mod);
        }
        return ret;
    }

};


#include "dmrg/models/chem/model_qc.hpp"

#endif
