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

#ifndef REL_QC_MODEL_H
#define REL_QC_MODEL_H

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

#include "dmrg/models/chem/rel/term_maker.h"
#include "dmrg/models/chem/rel/rel_chem_detail.h"
#include "dmrg/models/chem/pg_util.h"

template<class Matrix, class SymmGroup>
class rel_qc_model : public model_impl<Matrix, SymmGroup>
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
    
    rel_qc_model(Lattice const & lat_, BaseParameters & parms_);
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    // align function moved here from rel_chem_detail.h
    static rel_chem_detail::IndexTuple align(int i, int j, int k, int l) {
        // For the time being (12.06.14) do nothing in the relativistic model
        return rel_chem_detail::IndexTuple(i,j,k,l);
    }

    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        // type == site for lattice = spinors
        return phys_indices[type];
    }
    tag_type identity_matrix_tag(size_t type) const
    {
		// transfer operator
        //if (type == (lat.size()/2) -1)
        //    return ident_transfer;
        if (type < lat.size()/2)
            return ident_unbar;
        else
            return ident_bar;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
		// transfer operator
        //if (type == (lat.size()/2) -1)
        //    return fill_transfer;
        if (type < lat.size()/2)
            return fill_unbar;
        else
            return fill_bar;
    }

    bool is_term_allowed(int i, int j, int k, int l)
    {
        typename SymmGroup::charge I(0), J(0), K(0), L(0), tmp(0);
        typename SymmGroup::charge charges[] = {I,J,K,L};
        std::size_t site[] = {i, j, k, l};
        for (int ii=0; ii<4; ++ii) {
            charges[ii][2] = lat.get_prop<int>("irrep", site[ii]);
            if (site[ii] < lat.size()/2) {charges[ii][0] = 1;}
            else if (site[ii] >= lat.size()/2) {charges[ii][1] = 1;}
            else {throw std::runtime_error("integrals parsing failed\n");}
        if (ii%2 == 0) {
            tmp = SymmGroup::fuse(tmp,charges[ii]);}
        else if (ii%2 == 1) {
            tmp = SymmGroup::fuse(tmp, -charges[ii]);}
        //maquis::cout << "site: " << site[ii] << " charge: " << charges[ii] << std::endl;
        }
        //maquis::cout << "(" << i << j << k << l << "): " << tmp << std::endl;
        if (tmp[0] == 0 && tmp[1] == 0 &&  tmp[2] != 0) {return false;}
        else {return true;}
    }


    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return rel_chem_detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "create_unbar")
            return create_unbar;
        else if (name == "create_bar")
            return create_bar;
        else if (name == "destroy_unbar")
            return destroy_unbar;
        else if (name == "destroy_bar")
            return destroy_bar;
        else if (name == "count_unbar")
            return count_unbar;
        else if (name == "count_bar")
            return count_bar;
        //else if (name == "e2d")
        //    return e2d;
        //else if (name == "d2e")
        //    return d2e;
        //else if (name == "docc")
        //    return docc;
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

        op_t create_unbar_op, create_bar_op, destroy_unbar_op, destroy_bar_op,
             count_unbar_op, count_bar_op,
             //docc_op, e2d_op, d2e_op,
             //swap_b2u_op, swap_u2b_op,
             //create_unbar_count_bar_op, create_bar_count_unbar_op, destroy_unbar_count_bar_op, destroy_bar_count_unbar_op,
             ident_unbar_op, ident_bar_op, fill_unbar_op, fill_bar_op;

        ident_unbar_op = tag_handler->get_op(ident_unbar);
        ident_bar_op = tag_handler->get_op(ident_bar);
        fill_unbar_op = tag_handler->get_op(fill_unbar);
        fill_bar_op = tag_handler->get_op(fill_bar);
        create_unbar_op = tag_handler->get_op(create_unbar);
        create_bar_op = tag_handler->get_op(create_bar);
        destroy_unbar_op = tag_handler->get_op(destroy_unbar);
        destroy_bar_op = tag_handler->get_op(destroy_bar);
        count_unbar_op = tag_handler->get_op(count_unbar);
        count_bar_op = tag_handler->get_op(count_bar);
        //e2d_op = tag_handler->get_op(e2d);
        //d2e_op = tag_handler->get_op(d2e);
        //docc_op = tag_handler->get_op(docc);

        //gemm(destroy_bar_op, create_unbar_op, swap_b2u_op); // S_plus
        //gemm(destroy_unbar_op, create_bar_op, swap_u2b_op); // S_minus

        //gemm(count_bar_op, create_unbar_op, create_unbar_count_bar_op);
        //gemm(count_unbar_op, create_bar_op, create_bar_count_unbar_op);
        //gemm(count_bar_op, destroy_unbar_op, destroy_unbar_count_bar_op);
        //gemm(count_unbar_op, destroy_bar_op, destroy_bar_count_unbar_op);

        #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

        //GENERATE_SITE_SPECIFIC(ident_unbar_op)
        //GENERATE_SITE_SPECIFIC(ident_bar_op)
        //GENERATE_SITE_SPECIFIC(fill_unbar_op)
        //GENERATE_SITE_SPECIFIC(fill_bar_op)
        std::vector<op_t> ident_ops = this->generate_site_specific_ops(ident_unbar_op, ident_bar_op);
        std::vector<op_t> fill_ops = this->generate_site_specific_ops(fill_unbar_op, fill_bar_op);
        std::vector<op_t> create_ops = this->generate_site_specific_ops(create_unbar_op, create_bar_op);
        std::vector<op_t> destroy_ops = this->generate_site_specific_ops(destroy_unbar_op, destroy_bar_op);
        std::vector<op_t> count_ops = this->generate_site_specific_ops(count_unbar_op, count_bar_op);
        
//         GENERATE_SITE_SPECIFIC(create_unbar_op)
//         GENERATE_SITE_SPECIFIC(create_bar_op)
//         GENERATE_SITE_SPECIFIC(destroy_unbar_op)
//         GENERATE_SITE_SPECIFIC(destroy_bar_op)
//         GENERATE_SITE_SPECIFIC(count_unbar_op)
//         GENERATE_SITE_SPECIFIC(count_bar_op)
        std::vector<op_t> create_unbar_ops = this->generate_site_specific_ops(create_unbar_op, fill_bar_op);
		std::vector<op_t> create_bar_ops = this->generate_site_specific_ops(fill_unbar_op, create_bar_op);
		std::vector<op_t> destroy_unbar_ops = this->generate_site_specific_ops(destroy_unbar_op, fill_bar_op);
		std::vector<op_t> destroy_bar_ops = this->generate_site_specific_ops(fill_unbar_op, destroy_bar_op);
		std::vector<op_t> count_unbar_ops = this->generate_site_specific_ops(count_unbar_op, ident_bar_op);
		std::vector<op_t> count_bar_ops = this->generate_site_specific_ops(ident_unbar_op, count_bar_op);


        // All this cases make no sense using a spinor lattice!
        //GENERATE_SITE_SPECIFIC(e2d_op)
        //GENERATE_SITE_SPECIFIC(d2e_op)
        //GENERATE_SITE_SPECIFIC(docc_op)

        //GENERATE_SITE_SPECIFIC(swap_b2u_op)
        //GENERATE_SITE_SPECIFIC(swap_u2b_op)
        //GENERATE_SITE_SPECIFIC(create_unbar_count_bar_op)
        //GENERATE_SITE_SPECIFIC(create_bar_count_unbar_op)
        //GENERATE_SITE_SPECIFIC(destroy_unbar_count_bar_op)
        //GENERATE_SITE_SPECIFIC(destroy_bar_count_unbar_op)

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
                    if (it->value() == "Nunbar")
                        meas_op = count_unbar_ops;
                    else if (it->value() == "Nbar")
                        meas_op = count_bar_ops;
                    else if (it->value() == "N")
                        meas_op = count_ops;
                    //else if (it->value() == "Nunbar*Nbar" || it->value() == "docc")
                    //    meas_op = docc_ops;
                    else
                        throw std::runtime_error("Invalid observable\nLocal measurements supported so far are \"Nunbar\" and \"Nbar\"\n");

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

                std::string bra_ckp("");
                if(lhs == "MEASURE_TRANSITION_TWOPTDM"){
                    name = "transition_twoptdm";
                    bra_ckp = it->value();
                }
                else
                    name = "twoptdm";

                std::vector<bond_element> synchronous_meas_operators;
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(create_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_unbar_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(create_bar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_bar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_unbar_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_bar_ops, true) );
                meas_operators.push_back( std::make_pair(create_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_unbar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_bar_ops, true) );
                synchronous_meas_operators.push_back(meas_operators);
                }
                {
                bond_element meas_operators;
                meas_operators.push_back( std::make_pair(create_bar_ops, true) );
                meas_operators.push_back( std::make_pair(create_bar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_bar_ops, true) );
                meas_operators.push_back( std::make_pair(destroy_bar_ops, true) );
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
                    if (*it2 == "c_dag") {
                        meas_operators.push_back( std::make_pair(create_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c") {
                        meas_operators.push_back( std::make_pair(destroy_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "N") {
                        meas_operators.push_back( std::make_pair(count_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_unbar") {
                        meas_operators.push_back( std::make_pair(destroy_unbar_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "c_bar") {
                        meas_operators.push_back( std::make_pair(destroy_bar_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_unbar") {
                        meas_operators.push_back( std::make_pair(create_unbar_ops, true) );
                        ++f_ops;
                    }
                    else if (*it2 == "cdag_bar") {
                        meas_operators.push_back( std::make_pair(create_bar_ops, true) );
                        ++f_ops;
                    }

                    else if (*it2 == "id" || *it2 == "Id") {
                        meas_operators.push_back( std::make_pair(ident_ops, false) );
                    }
                    else if (*it2 == "Nunbar") {
                        meas_operators.push_back( std::make_pair(count_unbar_ops, false) );
                    }
                    else if (*it2 == "Nbar") {
                        meas_operators.push_back( std::make_pair(count_bar_ops, false) );
                    }
                    else if (*it2 == "docc" || *it2 == "Nunbar*Nbar") {
                    //    meas_operators.push_back( std::make_pair(docc_ops, false) );
                    }
                    else if (*it2 == "cdag_unbar*c_bar" || *it2 == "splus") {
                    //    meas_operators.push_back( std::make_pair(swap_b2u_ops, false) );
                    }
                    else if (*it2 == "cdag_bar*c_unbar" || *it2 == "sminus") {
                    //    meas_operators.push_back( std::make_pair(swap_u2b_ops, false) );
                    }

                    else if (*it2 == "cdag_unbar*cdag_bar") {
                    //    meas_operators.push_back( std::make_pair(e2d_ops, false) );
                    }
                    else if (*it2 == "c_unbar*c_bar") {
                    //    meas_operators.push_back( std::make_pair(d2e_ops, false) );
                    }

                    else if (*it2 == "cdag_unbar*Nbar") {
                    //    meas_operators.push_back( std::make_pair(create_unbar_count_bar_ops, true) );
                    //    ++f_ops;
                    }
                    else if (*it2 == "cdag_bar*Nunbar") {
                    //    meas_operators.push_back( std::make_pair(create_bar_count_unbar_ops, true) );
                    //    ++f_ops;
                    }
                    else if (*it2 == "c_unbar*Nbar") {
                    //    meas_operators.push_back( std::make_pair(destroy_unbar_count_bar_ops, true) );
                    //    ++f_ops;
                    }
                    else if (*it2 == "c_bar*Nunbar") {
                    //    meas_operators.push_back( std::make_pair(destroy_bar_count_unbar_ops, true) );
                    //    ++f_ops;
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
//     Index<SymmGroup> phys;
    std::vector<Index<SymmGroup> > phys_indices;

    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    tag_type ident_unbar, ident_bar, fill_unbar, fill_bar,
             create_unbar, create_bar, destroy_unbar, destroy_bar,
             count_unbar, count_bar, ident_transfer, fill_transfer;

    typename SymmGroup::subcharge max_irrep;

    std::vector<op_t> generate_site_specific_ops(op_t const & op_unbar, op_t const & op_bar) const
    {
        PGDecorator<SymmGroup> set_symm;
        std::vector<op_t> ret;
        for (std::size_t site = 0; site < lat.size()/2; ++site) {
            int irrep = lat.get_prop<int>("irrep", site);
            //for (typename SymmGroup::subcharge irrep = 0; irrep < max_irrep + 1; ++irrep) {
                op_t mod(set_symm(op_unbar.basis(), irrep));
                for (std::size_t b = 0; b < op_unbar.n_blocks(); ++b)
                    mod[b] = op_unbar[b];

                //maquis::cout << mod << std::endl;
                ret.push_back(mod);
            //}
        }
        for (std::size_t site = lat.size()/2; site < lat.size(); ++site) {
            int irrep = lat.get_prop<int>("irrep", site);
            //for (typename SymmGroup::subcharge irrep = 0; irrep < max_irrep + 1; ++irrep) {
                op_t mod(set_symm(op_bar.basis(), irrep));
                for (std::size_t b = 0; b < op_bar.n_blocks(); ++b)
                    mod[b] = op_bar[b];
                
                //maquis::cout << mod << std::endl;
                ret.push_back(mod);
            //}
        }
        return ret;
    }
};


#include "dmrg/models/chem/rel/rel_model_qc.hpp"

#endif
