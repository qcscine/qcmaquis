/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef APP_ALPS_MODEL_H
#define APP_ALPS_MODEL_H

#include "dmrg/models/model.h"

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include "symm_handler.hpp"


template <class Matrix, class SymmGroup>
class ALPSModel : public Model<Matrix, SymmGroup>
{
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
    typedef Hamiltonian<Matrix, SymmGroup> ham;
    typedef typename ham::table_type table_type;
    typedef typename ham::table_ptr table_ptr;
    typedef typename table_type::tag_type tag_type;
    typedef typename std::pair<value_type, tag_type> term_descriptor;
    
    typedef typename ham::op_t op_t;
    typedef typename ham::hamtagterm_t hamtagterm_t;
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;

    typedef std::pair<std::string, int> opkey_type;
    typedef std::map<opkey_type, tag_type> opmap_type;
    typedef typename opmap_type::const_iterator opmap_const_iterator;
    
    typedef typename SymmGroup::charge charge;
    
    
    ALPSModel (const graph_type& lattice_, const alps::Parameters& parms_)
    : parms(parms_)
    , lattice(lattice_)
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
        std::vector<std::vector<term_descriptor> > site_terms( alps::maximum_vertex_type(lattice.graph())+1 );
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

                hamtagterm_t term;
                term.scale = 1.;
                term.with_sign = false;
                term.fill_operator = operators[opkey_type("ident", type)];
                term.operators.push_back( std::make_pair(p, match->second) );
                terms.push_back(term);
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
                
                opmap_const_iterator match1 = operators.find(opkey_type(simplify_name(op1), type_s));
                if (match1 == operators.end())
                    match1 = register_operator(op1, type_s, parms);
                opmap_const_iterator match2 = operators.find(opkey_type(simplify_name(op2), type_t));
                if (match2 == operators.end())
                    match2 = register_operator(op2, type_t, parms);

                bool with_sign = fermionic(b1, op1, b2, op2);
                
                hamtagterm_t term;
                term.scale = boost::get<0>(*tit).value();
                term.with_sign = with_sign;
                term.fill_operator = operators[opkey_type( (with_sign) ? "fill" : "ident", type_s)];
                
                {
                    tag_type mytag = match1->second;
                    if (with_sign && !wrap_pbc) {
                        // Note inverse notation because of notation in operator.
                        std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_s)],
                                                                                            mytag);
                        mytag = ptag.first;
                        term.scale *= ptag.second;
                    }
                    if (with_sign && wrap_pbc)
                        term.scale *= -1.;
                    term.operators.push_back( std::make_pair(p_s, mytag) );
                }
                {
                    tag_type mytag = match2->second;
                    if (with_sign && wrap_pbc) {
                        // Note inverse notation because of notation in operator.
                        std::pair<tag_type, value_type> ptag = tag_handler->get_product_tag(operators[opkey_type("fill",type_t)],
                                                                                            mytag);
                        mytag = ptag.first;
                        term.scale *= ptag.second;
                    }
                    term.operators.push_back( std::make_pair(p_t, mytag) );
                }
                
                terms.push_back(term);
            }
        }
    }
            
    op_t const & get_identity() const
    {
        return tag_handler->get_op( operators[opkey_type("ident", 0)] );
    }
    
    Index<SymmGroup> get_phys() const
    {
        return symm_basis[0].phys();
    }
    
    Hamiltonian<Matrix, SymmGroup> H() const
    {
        std::vector<typename ham::hamterm_t> terms_ops;
        return ham(get_phys(), get_identity(), terms_ops, operators[opkey_type("ident", 0)], terms, tag_handler);
    }
    
    typename SymmGroup::charge initc (BaseParameters& parms_) const
    {
        return init_charge<SymmGroup>(parms_, all_conserved_qn);
    }
    
    block_matrix<Matrix, SymmGroup> get_op(std::string const & name) const
    {
        // TODO: hard coding site type = 0 !!
        int type = 0;

        if (name == "id" || name == "ident" || name == "identity") {
            return tag_handler->get_op( operators[opkey_type("ident", type)] );
        } else {
            opmap_const_iterator match = operators.find(opkey_type(name, type));
            if (match == operators.end()) {
                SiteOperator op = make_site_term(name, parms);
                match = register_operator(op, type, parms);
            }
            return tag_handler->get_op(match->second);
        }
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
    
    alps::Parameters const & parms;
    const graph_type& lattice;
    alps::model_helper<I> model;
    mutable table_ptr tag_handler;

    std::vector<symmetric_basis_descriptor<SymmGroup> > symm_basis;
    
    mutable opmap_type operators; // key=<name,type>

    std::vector<hamtagterm_t> terms;
    qn_map_type all_conserved_qn;
};

// Loading Measurements
template <class Matrix, class SymmGroup>
Measurements<Matrix, SymmGroup> ALPSModel<Matrix, SymmGroup>::measurements () const
{
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    Measurements<Matrix, SymmGroup> meas;
    
    // TODO: hard coding site type = 0 !!
    int type = 0;
    
    meas.set_identity(get_op("ident"));
    
    {
        boost::regex average_expr("^MEASURE_AVERAGE\\[(.*)]$");
        boost::regex locale_expr("^MEASURE_LOCAL\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            mterm_t term;
            
            bool found = false;
            if (boost::regex_match(lhs, what, average_expr)) {
                term.type = mterm_t::Average;
                found = true;
            }
            if (boost::regex_match(lhs, what, locale_expr)) {
                term.type = mterm_t::Local;
                found = true;
            }
            
            if (found) {
				alps::SiteBasisDescriptor<I> b = model.site_basis(type);

                if (model.has_bond_operator(it->value())) {
                	BondOperator bondop = model.get_bond_operator(it->value());

                    typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
                    V  ops = bondop.template templated_split<value_type>(b,b);
                    for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                        SiteOperator op1 = tit->template get<1>();
                        SiteOperator op2 = tit->template get<2>();

                        bool with_sign = fermionic(b, op1, b, op2);

                        std::ostringstream ss;
                        ss << what.str(1);
                        if (ops.size() > 1) ss << " (" << int(tit-ops.begin())+1 << ")";
                        term.name = ss.str();

                        term.fill_operator = (with_sign) ? get_op("fill") : get_op("ident");
                        {
                        	op_t tmp;
                        	if (with_sign)
                        		gemm(term.fill_operator, get_op(simplify_name(op1)), tmp); // Note inverse notation because of notation in operator.
                            else
                                tmp = get_op(simplify_name(op1));
                        	term.operators.push_back( std::make_pair(tit->template get<0>().value()*tmp, b.is_fermionic(simplify_name(op1))) );
                        }
                        {
                            term.operators.push_back( std::make_pair(get_op(simplify_name(op2)), b.is_fermionic(simplify_name(op2))) );
                        }
	                    meas.add_term(term);
                    }
                } else {
                    term.name = what.str(1);

					SiteOperator op = make_site_term(it->value(), parms);

					if (b.is_fermionic(simplify_name(op)))
						throw std::runtime_error("Cannot measure local fermionic operators.");

					term.operators.push_back( std::make_pair(get_op(it->value()), false) );
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
                alps::SiteBasisDescriptor<I> b = model.site_basis(type);
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
                

                boost::char_separator<char> op_sep(":");
                tokenizer ops_tokens(parts[0], op_sep);
                for (tokenizer::iterator it_op=ops_tokens.begin();
                     it_op != ops_tokens.end(); it_op++)
                {
                    SiteOperator op = make_site_term(*it_op, parms);
                    bool f = b.is_fermionic(simplify_name(op));
                    term.operators.push_back( std::make_pair(get_op(simplify_name(op)), f) );
                    if (f) ++f_ops;
                }

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
                                
                term.fill_operator = (f_ops > 0) ? get_op("fill") : get_op("ident");
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
            
            mterm_t term;
            term.fill_operator = get_op("ident");
            
            if (boost::regex_match(lhs, what, expression)) {
                value = it->value();
                term.name = what.str(1);
                term.type = mterm_t::MPSBonds;
                
                alps::SiteBasisDescriptor<I> b = model.site_basis(type);
                int f_ops = 0;
                
                boost::char_separator<char> sep(":");
                tokenizer corr_tokens(value, sep);
                for (tokenizer::iterator it2=corr_tokens.begin();
                     it2 != corr_tokens.end();
                     it2++)
                {
                    SiteOperator op = make_site_term(*it2, parms);
                    bool f = b.is_fermionic(simplify_name(op));
                    term.operators.push_back( std::make_pair(get_op(simplify_name(op)), f) );
                    if (f) ++f_ops;
                }
                if (term.operators.size() == 1) {
                    term.operators.push_back(term.operators[0]);
                    if (term.operators[1].second) ++f_ops;
                }
                
                if (f_ops > 0) {
                    term.fill_operator = get_op("fill");
                }
                

                if (f_ops % 2 != 0)
                    throw std::runtime_error("Number of fermionic operators has to be even.");

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
            term.fill_operator = get_op("ident");
            
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
                
                alps::SiteBasisDescriptor<I> b = model.site_basis(type);
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
                    SiteOperator op = make_site_term(*it2, parms);
                    bool f = b.is_fermionic(simplify_name(op));
                    term.operators.push_back( std::make_pair(get_op(simplify_name(op)), f) );
                    if (f) ++f_ops;
                }
                if (term.operators.size() == 1) {
                    term.operators.push_back(term.operators[0]);
                    if (term.operators[1].second) ++f_ops;
                }
                
                
                if (f_ops > 0) {
                    term.fill_operator = get_op("fill");
                }

                if (f_ops % 2 != 0)
                    throw std::runtime_error("Number of fermionic operators has to be even.");
                
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
