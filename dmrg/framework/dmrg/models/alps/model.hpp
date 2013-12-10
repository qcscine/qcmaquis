/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef APP_ALPS_MODEL_H
#define APP_ALPS_MODEL_H

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/models/alps/lattice.hpp"

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include "symm_handler.hpp"

namespace detail {
    inline alps::graph_helper<> const& get_graph(Lattice const& lat_)
    {
        alps_lattice const* alattice = static_cast<alps_lattice const*>(lat_.impl().get());
        return alattice->alps_graph();
    }
}

template <class Matrix, class SymmGroup>
class ALPSModel : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;
    
    typedef alps::SiteOperator SiteOperator;
    typedef alps::BondOperator BondOperator;
    
    typedef typename Matrix::value_type value_type;
    typedef typename maquis::traits::scalar_type<Matrix>::type scalar_type;
    typedef boost::multi_array<value_type,2> alps_matrix;
    typedef std::map<std::string, int> qn_map_type;
    
    typedef short I;
    typedef alps::graph_helper<> graph_type;
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    
public:
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
    typedef typename base::size_t size_t;
    
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;

    typedef std::pair<std::string, int> opkey_type;
    typedef std::map<opkey_type, tag_type> opmap_type;
    typedef typename opmap_type::const_iterator opmap_const_iterator;
    
    typedef typename SymmGroup::charge charge;
    
    
    ALPSModel (Lattice const& lattice_, const alps::Parameters& parms_)
    : parms(parms_)
    , lattice(detail::get_graph(lattice_))
    , model(lattice, parms)
    , tag_handler(new table_type())
    {
        locale_shared i;
        
        symm_basis.resize( alps::maximum_vertex_type(lattice.graph())+1 );
        
        /// Parsing conserved quantum numbers
        std::set<std::string> all_qn;
        for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
            std::set<std::string> type_qn = model.quantum_numbers(type);
            all_qn.insert(type_qn.begin(), type_qn.end());
        }

        if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
            boost::char_separator<char> sep(" ,");
            std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
            tokenizer qn_tokens(qn_string, sep);
            int n=0;
            for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
                if (parms.defined(*it + "_total")) {
                    if (all_qn.find(*it) != all_qn.end())
                        all_conserved_qn.insert( std::make_pair(*it, n++) );
                    else
                        throw std::runtime_error("quantumnumber "+(*it)+" not defined in the model.");
                }
            }
        }
        
        /// Load all possible basis
        for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
            alps::SiteBasisDescriptor<I> b = model.site_basis(type);
            alps::site_basis<I> states(b);
            symm_basis[type] = symmetric_basis_descriptor<SymmGroup>(b, all_conserved_qn);
            
            op_t ident, fill;
            for (int i=0; i<symm_basis[type].size(); ++i) {
                charge c = symm_basis[type].charge(i);
                size_t bsize = symm_basis[type].block_size(i);
                // maquis::cout << "Inserting " << c << " for " << states[i] << std::endl;
                
                if (!ident.has_block(c, c))
                    ident.insert_block(Matrix::identity_matrix(bsize), c, c);
                
                int sign = (alps::is_fermionic(b, states[i])) ? -1 : 1;
                if (!fill.has_block(c, c))
                    fill.insert_block(Matrix::identity_matrix(bsize), c, c);
                fill(symm_basis[type].coords(i), symm_basis[type].coords(i)) = sign;
            }
            operators[opkey_type("ident", type)] = tag_handler->register_op(ident, tag_detail::bosonic);
            operators[opkey_type("fill",  type)] = tag_handler->register_op(fill,  tag_detail::bosonic);
        }
        
        
        /// site_term loop with cache to avoid recomputing matrices
        std::vector<std::vector<std::pair<value_type, tag_type> > > site_terms( alps::maximum_vertex_type(lattice.graph())+1 );
        for (graph_type::site_iterator it=lattice.sites().first; it!=lattice.sites().second; ++it) {
            int p = lattice.vertex_index(*it);
            int type = lattice.site_type(*it);
            
            if (site_terms[type].size() == 0) {
                typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator> > V;
                V  ops = model.site_term(type).template templated_split<value_type>();
                                                        
                for (int n=0; n<ops.size(); ++n) {
                    if (boost::get<0>(ops[n]).value() != 0.) {
                        SiteOperator op = boost::get<1>(ops[n]);
                        opmap_const_iterator match = operators.find(opkey_type(simplify_name(op), type));
                        if (match == operators.end())
                            match = register_operator(op, type, parms);
                        site_terms[type].push_back( std::make_pair(boost::get<0>(ops[n]).value(), match->second)  );
                    }
                }
            }

            // All site terms summed into one
            if (site_terms[type].size() > 0) {
                opmap_const_iterator match = operators.find(opkey_type("site_terms", type));
                if (match == operators.end()) {
                    op_t op_matrix;
                    for (int n=0; n<site_terms[type].size(); ++n)
                        op_matrix += site_terms[type][n].first * tag_handler->get_op(site_terms[type][n].second);
                    tag_type mytag = tag_handler->register_op(op_matrix, tag_detail::bosonic);
                    boost::tie(match, boost::tuples::ignore) = operators.insert( std::make_pair(opkey_type("site_terms", type), mytag) );
                }

                term_descriptor term;
                term.coeff = 1.;
                term.is_fermionic = false;
                term.push_back( boost::make_tuple(p, match->second) );
                this->terms_.push_back(term);
            }
        }
        
        /// bond terms loop
        for (graph_type::bond_iterator it=lattice.bonds().first; it!=lattice.bonds().second; ++it) {
            int p_s = lattice.source(*it);
            int p_t = lattice.target(*it);
            int type = lattice.bond_type(*it);
            int type_s = lattice.site_type(lattice.source(*it));
            int type_t = lattice.site_type(lattice.target(*it));
            
            bool wrap_pbc = boost::get(alps::boundary_crossing_t(), lattice.graph(), *it);
            
            BondOperator bondop = model.bond_term(type);
            
            typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
            alps::SiteBasisDescriptor<I> b1 = model.site_basis(type_s);
            alps::SiteBasisDescriptor<I> b2 = model.site_basis(type_t);
            
            
            V  ops = bondop.template templated_split<value_type>(b1,b2);
            for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                SiteOperator op1 = boost::get<1>(*tit);
                SiteOperator op2 = boost::get<2>(*tit);
                
                if (alps::numeric::is_nonzero(boost::get<0>(*tit).value())) {
                    
                    opmap_const_iterator match1 = operators.find(opkey_type(simplify_name(op1), type_s));
                    if (match1 == operators.end())
                        match1 = register_operator(op1, type_s, parms);
                    opmap_const_iterator match2 = operators.find(opkey_type(simplify_name(op2), type_t));
                    if (match2 == operators.end())
                        match2 = register_operator(op2, type_t, parms);
                    
                    bool with_sign = fermionic(b1, op1, b2, op2);
                    
                    term_descriptor term;
                    term.coeff = boost::get<0>(*tit).value();
                    term.is_fermionic = with_sign;
                    
                    {
                        tag_type mytag = match1->second;
                        if (with_sign && !wrap_pbc) {
                            // Note inverse notation because of notation in operator.
                            std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_s)],
                                                                                                mytag);
                            mytag = ptag.first;
                            term.coeff *= ptag.second;
                        }
                        if (with_sign && wrap_pbc)
                            term.coeff *= -1.;
                        term.push_back( boost::make_tuple(p_s, mytag) );
                    }
                    {
                        tag_type mytag = match2->second;
                        if (with_sign && wrap_pbc) {
                            // Note inverse notation because of notation in operator.
                            std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_t)],
                                                                                                mytag);
                            mytag = ptag.first;
                            term.coeff *= ptag.second;
                        }
                        term.push_back( boost::make_tuple(p_t, mytag) );
                    }
                    
                    this->terms_.push_back(term);
                }
            }
        }
    }
    
    Index<SymmGroup> const& phys_dim(size_t type) const
    {
        return symm_basis[type].phys_dim();
    }
    
    typename SymmGroup::charge total_quantum_numbers(BaseParameters& parms_) const
    {
        return init_charge<SymmGroup>(parms_, all_conserved_qn);
    }
    
    tag_type identity_matrix_tag(size_t type) const
    {
        return operators[opkey_type("ident", type)]; // TODO: avoid using map here
    }

    tag_type filling_matrix_tag(size_t type) const
    {
        return operators[opkey_type("fill", type)]; // TODO: avoid using map here
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "id" || name == "ident" || name == "identity") {
            return operators[opkey_type("ident", type)];
        } else {
            opmap_const_iterator match = operators.find(opkey_type(name, type));
            if (match == operators.end()) {
                SiteOperator op = make_site_term(name, parms);
                match = register_operator(op, type, parms);
            }
            return match->second;
        }
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }

    Measurements<Matrix, SymmGroup> measurements () const;

private:
    
    template <class SiteOp>
    std::string simplify_name(const SiteOp &op) const
    {
        std::string term = op.term();
        std::string arg = "("+op.site()+")";
        boost::algorithm::replace_all(term,arg,"");
        return term;
    }
    
    bool fermionic (alps::SiteBasisDescriptor<I> b1, SiteOperator op1,
                    alps::SiteBasisDescriptor<I> b2, SiteOperator op2) const
    {
        return b1.is_fermionic(simplify_name(op1)) && b2.is_fermionic(simplify_name(op2));
    }
    
    inline op_t convert_matrix (const alps_matrix& m, int type) const
    {
        op_t newm;
        for (int i=0; i<m.shape()[0]; ++i) {
            for (int j=0; j<m.shape()[1]; ++j) {
                if (m[i][j] != 0.) {
                    charge c_i = symm_basis[type].charge(i);
                    size_t bsize_i = symm_basis[type].block_size(i);
                    charge c_j = symm_basis[type].charge(j);
                    size_t bsize_j = symm_basis[type].block_size(j);
                    
                    if (!newm.has_block(c_i, c_j))
                        newm.insert_block(Matrix(bsize_i, bsize_j, 0), c_i, c_j);
                    // Notation: going from state i to state j
                    newm(symm_basis[type].coords(i), symm_basis[type].coords(j)) = m[i][j];
                }
            }
        }
        return newm;
    }
    
    alps::SiteOperator make_site_term(std::string x, alps::Parameters const & parms) const
    {
        if (x[x.size()-1]!=')')
            x += "(i)";
        alps::SiteOperator op(x,"i");
        model.substitute_operators(op, parms);
        return op;
    }
    
    opmap_const_iterator register_operator(SiteOperator const& op, int type, alps::Parameters const& p) const
    {
        alps::SiteBasisDescriptor<I> b = model.site_basis(type);
        alps_matrix m = alps::get_matrix(value_type(), op, b, p, true);
        tag_detail::operator_kind kind = b.is_fermionic(simplify_name(op)) ? tag_detail::fermionic : tag_detail::bosonic;
        tag_type mytag = tag_handler->register_op(convert_matrix(m, type), kind);
        
        opmap_const_iterator match;
        boost::tie(match, boost::tuples::ignore) = operators.insert( std::make_pair(opkey_type(simplify_name(op), type), mytag) );
        return match;
    }
    
    void fill_meas_with_operators(mterm_t & term, std::string const& ops, bool repeat_one=false) const
    {
        int ntypes = alps::maximum_vertex_type(lattice.graph())+1;
        int f_ops = 0;
        
        boost::char_separator<char> sep(":");
        tokenizer corr_tokens(ops, sep);
        for (tokenizer::iterator it2=corr_tokens.begin(); it2 != corr_tokens.end(); it2++)
        {
            enum {uknown, bosonic, fermionic} kind = uknown;
            std::vector<op_t> tops(ntypes);
            for (int type=0; type<ntypes; ++type) {
                alps::SiteBasisDescriptor<I> b = model.site_basis(type);
                if (b.has_operator(*it2)) {
                    SiteOperator op = make_site_term(*it2, parms);
                    bool is_ferm = b.is_fermionic(simplify_name(op));
                    if (kind == uknown)
                        kind = is_ferm ? fermionic : bosonic;
                    else if ((is_ferm && kind==bosonic) || (!is_ferm && kind==fermionic))
                        throw std::runtime_error("Model is inconsitent. On some site the operator " + *it2 + "fermionic, on others is bosonic.");
                    tops[type] = this->get_operator(*it2, type);
                }
            }
            
            if (kind == fermionic) ++f_ops;
            term.operators.push_back( std::make_pair(tops, (kind==fermionic)) );
        }
        if (repeat_one && term.operators.size() == 1) {
            term.operators.push_back(term.operators[0]);
            if (term.operators[1].second) ++f_ops;
        }
        
        if (f_ops % 2 != 0)
            throw std::runtime_error("Number of fermionic operators has to be even.");
    }
    
    alps::Parameters const& parms;
    graph_type const& lattice;
    alps::model_helper<I> model;
    mutable table_ptr tag_handler;

    std::vector<symmetric_basis_descriptor<SymmGroup> > symm_basis;
    
    mutable opmap_type operators; // key=<name,type>

    qn_map_type all_conserved_qn;
};

// Loading Measurements
template <class Matrix, class SymmGroup>
Measurements<Matrix, SymmGroup> ALPSModel<Matrix, SymmGroup>::measurements () const
{
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    
    int ntypes = alps::maximum_vertex_type(lattice.graph())+1;
    
    std::vector<op_t> identitities(ntypes), fillings(ntypes);
    for (int type=0; type<ntypes; ++type) {
        identitities[type] = this->identity_matrix(type);
        fillings[type]     = this->filling_matrix(type);
    }
    Measurements<Matrix, SymmGroup> meas(identitities, fillings);
    
    {
        boost::regex average_expr("^MEASURE_AVERAGE\\[(.*)]$");
        boost::regex locale_expr("^MEASURE_LOCAL\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            mterm_t term;
            std::string obsname;
            
            bool found = false;
            if (boost::regex_match(lhs, what, average_expr)) {
                term.type = mterm_t::Average;
                found = true;
                obsname = what.str(1);
            }
            if (boost::regex_match(lhs, what, locale_expr)) {
                term.type = mterm_t::Local;
                found = true;
                obsname = what.str(1);
            }
            
            if (found) {

                if (model.has_bond_operator(it->value())) {
                    throw std::runtime_error("not yet implemented!");
//                	BondOperator bondop = model.get_bond_operator(it->value());
//
//                    typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
//                    
//                    std::vector<op_t> tops(ntypes);
//                    for (int type=0; type<ntypes; ++type) {
//
//                    alps::SiteBasisDescriptor<I> b1 = model.site_basis(type1);
//                    alps::SiteBasisDescriptor<I> b2 = model.site_basis(type2);
//                    
//                    V  ops = bondop.template templated_split<value_type>(b,b);
//                    for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
//                        SiteOperator op1 = boost::get<1>(*tit);
//                        SiteOperator op2 = boost::get<2>(*tit);
//
//                        bool with_sign = fermionic(b1, op1, b2, op2);
//
//                        std::ostringstream ss;
//                        ss << what.str(1);
//                        if (ops.size() > 1) ss << " (" << int(tit-ops.begin())+1 << ")";
//                        term.name = ss.str();
//
//                        term.with_sign = with_sign;
//                        {
//                        	op_t tmp;
//                        	if (with_sign)
//                        		gemm(term.fill_operator, this->get_operator(simplify_name(op1), type), tmp); // Note inverse notation because of notation in operator.
//                            else
//                                tmp = this->get_operator(simplify_name(op1), type);
//                        	term.operators.push_back( std::make_pair(boost::get<0>(*tit).value()*tmp, b.is_fermionic(simplify_name(op1))) );
//                        }
//                        {
//                            term.operators.push_back( std::make_pair(this->get_operator(simplify_name(op2), type), b.is_fermionic(simplify_name(op2))) );
//                        }
//	                    meas.add_term(term);
//                    }
                } else {
                    term.name = obsname;
                    
                    std::vector<op_t> tops(ntypes);
                    for (int type=0; type<ntypes; ++type) {
                        alps::SiteBasisDescriptor<I> b = model.site_basis(type);
                        if (b.has_operator(it->value())) {
                            SiteOperator op = make_site_term(it->value(), parms);
                            if (b.is_fermionic(simplify_name(op)))
                                throw std::runtime_error("Cannot measure local fermionic operators.");
                            
                            tops[type] = this->get_operator(it->value(), type);
                        }
                    }
					term.operators.push_back( std::make_pair(tops, false) );
                    meas.add_term(term);
                }
            }
        }
    }
    
    { // Example: MEASURE_LOCAL_AT[Custom correlation] = "bdag:b|(1,2),(3,4),(5,6)"
        boost::regex expression("^MEASURE_LOCAL_AT\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            std::string value = it->value();
            if (boost::regex_match(lhs, what, expression)) {
                int f_ops = 0;
                
                mterm_t term;
                term.type = mterm_t::LocalAt;
                term.name = what.str(1);
                
                boost::char_separator<char> part_sep("|");
                tokenizer part_tokens(value, part_sep);
                std::vector<std::string> parts;
                std::copy( part_tokens.begin(), part_tokens.end(), std::back_inserter(parts) );
                
                if (parts.size() != 2)
                    throw std::runtime_error("MEASURE_LOCAL_AT must contain a `|` delimiter.");
                
                /// parse operators
                fill_meas_with_operators(term, parts[0], false);

                /// parse positions
                boost::regex pos_re("\\(([^(^)]*)\\)");
                boost::sregex_token_iterator it_pos(parts[1].begin(), parts[1].end(), pos_re, 1);
                boost::sregex_token_iterator it_pos_end;
                for (; it_pos != it_pos_end; ++it_pos)
                {
                    boost::char_separator<char> int_sep(", ");
                    std::string raw = *it_pos;
                    tokenizer int_tokens(raw, int_sep);
                    
                    std::vector<std::size_t> pos;
                    std::transform(int_tokens.begin(), int_tokens.end(), std::back_inserter(pos),
                                   static_cast<std::size_t (*)(std::string const&)>(boost::lexical_cast<std::size_t, std::string>));
                    term.positions.push_back(pos);
                }
                if (f_ops % 2 != 0)
                    throw std::runtime_error("Number of fermionic operators has to be even.");
                
                meas.add_term(term);
            }
        }
    }

    {
        boost::regex expression("^MEASURE_MPS_BONDS\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            std::string value;
            
            if (boost::regex_match(lhs, what, expression)) {
                mterm_t term;
                term.name = what.str(1);
                term.type = mterm_t::MPSBonds;
                
                fill_meas_with_operators(term, it->value(), true);
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
                /// split op1:op2:...@p1,p2,p3,... into {op1:op2:...}, {p1,p2,p3,...}
                std::vector<std::string> value_split;
                boost::split( value_split, value, boost::is_any_of("@"));
                
                fill_meas_with_operators(term, value_split[0], true);

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

#endif
