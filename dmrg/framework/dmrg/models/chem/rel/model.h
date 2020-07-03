/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2015 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2015-2015 by Stefano Battaglia <stefabat@ethz.ch>
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
#include <regex>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/util.h"
#include "dmrg/models/chem/parse_integrals.h"
#include "dmrg/models/chem/pg_util.h"
#include "dmrg/models/chem/2u1/term_maker.h"
#include "dmrg/models/chem/rel/chem_helper.h"

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

    void create_terms();

    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }

    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        // type == site for lattice = spinors
        return phys_indices[type];
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident[type];
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill[type];
    }

    bool is_term_allowed(int i, int j, int k, int l)
    {
        typename SymmGroup::charge I(0), J(0), K(0), L(0), tmp(0);
        typename SymmGroup::charge charges[] = {I,J,K,L};
        std::size_t site[] = {i, j, k, l};
        for (int ii=0; ii<4; ++ii) {
            charges[ii][1] = lat.get_prop<int>("type", site[ii]);
            charges[ii][0] = 1;
        	if (ii%2 == 0) {
            	tmp = SymmGroup::fuse(tmp, charges[ii]);}
        	else if (ii%2 == 1) {
            	tmp = SymmGroup::fuse(tmp, -charges[ii]);}
        }

        if (tmp[0] == 0 && tmp[1] != parms["type"]) {return false;}
        else {return true;}
    }

    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return chem::detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "create")
            return create[type];
        else if (name == "destroy")
            return destroy[type];
        else if (name == "count")
            return count[type];
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

        op_t create_op, destroy_op, count_op,
             ident_op, fill_op;

        ident_op = tag_handler->get_op(ident[0]);
        fill_op = tag_handler->get_op(fill[0]);
        create_op = tag_handler->get_op(create[0]);
        destroy_op = tag_handler->get_op(destroy[0]);
        count_op = tag_handler->get_op(count[0]);

        #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

        GENERATE_SITE_SPECIFIC(ident_op)
        GENERATE_SITE_SPECIFIC(fill_op)
        GENERATE_SITE_SPECIFIC(create_op)
        GENERATE_SITE_SPECIFIC(destroy_op)
        GENERATE_SITE_SPECIFIC(count_op)

        #undef GENERATE_SITE_SPECIFIC

        measurements_type meas;

        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        {
            std::regex expression("^MEASURE_LOCAL\\[(.*)]$");
            std::smatch what;
            for (auto&& it: parms.get_range()) {
                std::string lhs = it.first;
                if (std::regex_match(lhs, what, expression)) {

                    op_vec meas_op;
                    if (it.second == "N")
                        meas_op = count_ops;
                    else
                        throw std::runtime_error("Invalid observable\nLocal measurement supported so far is \"N\"\n");

                    meas.push_back( new measurements::local<Matrix, SymmGroup>(what.str(1), lat, ident_ops, fill_ops, meas_op) );
                }
            }
        }

        {
        std::regex expression("^MEASURE_CORRELATIONS\\[(.*)]$");
        std::regex expression_half("^MEASURE_HALF_CORRELATIONS\\[(.*)]$");
        std::regex expression_nn("^MEASURE_NN_CORRELATIONS\\[(.*)]$");
        std::regex expression_halfnn("^MEASURE_HALF_NN_CORRELATIONS\\[(.*)]$");
        std::regex expression_oneptdm("^MEASURE\\[1rdm\\]");
        std::regex expression_transition_oneptdm("^MEASURE\\[trans1rdm\\]");
        std::regex expression_twoptdm("^MEASURE\\[2rdm\\]");
        std::regex expression_transition_twoptdm("^MEASURE\\[trans2rdm\\]");
        std::regex expression_threeptdm("^MEASURE\\[3rdm\\]");
        std::regex expression_transition_threeptdm("^MEASURE\\[trans3rdm\\]");
        std::smatch what;
        for (auto&& it : parms.get_range()) {
            std::string lhs = it.first;

            std::string name, value;
            bool half_only, nearest_neighbors_only;
            if (std::regex_match(lhs, what, expression)) {
                value = it.second;
                name = what.str(1);
                half_only = false;
                nearest_neighbors_only = false;
            }
            if (std::regex_match(lhs, what, expression_half)) {
                value = it.second;
                name = what.str(1);
                half_only = true;
                nearest_neighbors_only = false;
            }
            if (std::regex_match(lhs, what, expression_nn)) {
                value = it.second;
                name = what.str(1);
                half_only = false;
                nearest_neighbors_only = true;
            }
            if (std::regex_match(lhs, what, expression_halfnn)) {
                value = it.second;
                name = what.str(1);
                half_only = true;
                nearest_neighbors_only = true;
            }
            if (std::regex_match(lhs, what, expression_oneptdm) ||
                    std::regex_match(lhs, what, expression_transition_oneptdm)) {

                typedef std::vector<tag_type> tag_vec;
                typedef std::vector<tag_vec> bond_tag_element;
                typedef std::pair<std::vector<tag_vec>, value_type> scaled_bond_element;

                std::string bra_ckp("");
                if(lhs == "MEASURE[trans1rdm]"){
                    name = "transition_oneptdm";
                    bra_ckp = it.second;
                }
                else
                    name = "oneptdm";

	        	std::vector<scaled_bond_element> synchronous_meas_operators;
                {
                    bond_tag_element meas_operators;
                    meas_operators.push_back(create);
                    meas_operators.push_back(destroy);
                    synchronous_meas_operators.push_back(std::make_pair(meas_operators, 1));
                }
                half_only = true;
                std::vector<pos_t> positions;
                meas.push_back( new measurements::TaggedNRankRDM<Matrix, SymmGroup>(name, lat, tag_handler, ident, fill,
		                                                                            synchronous_meas_operators,half_only,
										                                            positions, bra_ckp));
            }

            else if (std::regex_match(lhs, what, expression_twoptdm) ||
                    std::regex_match(lhs, what, expression_transition_twoptdm)) {

                typedef std::vector<tag_type> tag_vec;
                typedef std::vector<tag_vec> bond_tag_element;
                typedef std::pair<std::vector<tag_vec>, value_type> scaled_bond_element;

                std::string bra_ckp("");
                if(lhs == "MEASURE[trans2rdm]"){
                    name = "transition_twoptdm";
                    bra_ckp = it.second;
                }
                else
                    name = "twoptdm";

		        std::vector<scaled_bond_element> synchronous_meas_operators;
                {
                    bond_tag_element meas_operators;
                    meas_operators.push_back(create);
                    meas_operators.push_back(create);
                    meas_operators.push_back(destroy);
                    meas_operators.push_back(destroy);
                    synchronous_meas_operators.push_back(std::make_pair(meas_operators, 1));
                }
                half_only = true;
                std::vector<pos_t> positions;
                meas.push_back( new measurements::TaggedNRankRDM<Matrix, SymmGroup>(name, lat, tag_handler, ident, fill,
		                                                                            synchronous_meas_operators,half_only,
										                                            positions, bra_ckp));
            }

            else if (!name.empty()) {
                typedef std::vector<tag_type> tag_vec;
                typedef std::vector<tag_vec> bond_tag_element;
                typedef std::pair<std::vector<tag_vec>, value_type> scaled_bond_element;

                int f_ops = 0;
                bond_tag_element meas_operators;

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
                        meas_operators.push_back(create);
                        ++f_ops;
                    }
                    else if (*it2 == "c") {
                        meas_operators.push_back(destroy);
                        ++f_ops;
                    }
                    else if (*it2 == "N") {
                        meas_operators.push_back(count);
                        ++f_ops;
                    }
                    else if (*it2 == "id" || *it2 == "Id") {
                        meas_operators.push_back(ident);
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

                std::vector<scaled_bond_element> synchronous_meas_operators;
                synchronous_meas_operators.push_back(std::make_pair(meas_operators, 1));
                meas.push_back( new measurements::TaggedNRankRDM<Matrix, SymmGroup>(name, lat, tag_handler, ident, fill, synchronous_meas_operators,
                                                                                    half_only, positions));
            }
        }
        }
        return meas;
    }

private:
    Lattice const & lat;
    BaseParameters & parms;
    std::vector<Index<SymmGroup> > phys_indices;

    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    std::vector<tag_type> ident, fill, create, destroy, count;

	std::vector<op_t> generate_site_specific_ops(op_t const & op) const
    {
        PGDecorator<SymmGroup> set_symm;
        std::vector<op_t> ret;
        for (typename SymmGroup::subcharge sc=0; sc < SymmGroup::get_max_irrep()+1; ++sc) {
            op_t mod(set_symm(op.basis(), sc));
            for (std::size_t b = 0; b < op.n_blocks(); ++b)
                mod[b] = op[b];

            ret.push_back(mod);
        }
        return ret;
    }

    std::vector<tag_type> register_site_specific(std::vector<op_t> const & ops, tag_detail::operator_kind kind)
    {
        std::vector<tag_type> ret;
        for (typename SymmGroup::subcharge sc=0; sc < SymmGroup::get_max_irrep()+1; ++sc) {
            std::pair<tag_type, value_type> newtag = tag_handler->checked_register(ops[sc], kind);
            assert( newtag.first < tag_handler->size() );
            assert( std::abs(newtag.second - value_type(1.)) == value_type() );
            ret.push_back(newtag.first);
        }

        return ret;
    }


};

#include "dmrg/models/chem/rel/model.hpp"

#endif
