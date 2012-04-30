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

template <class SymmGroup>
typename SymmGroup::charge convert_alps (alps::site_state<short> const & state, std::vector<std::pair<std::size_t, std::string> > const& qn);

template <class SymmGroup>
typename SymmGroup::charge init_charge (const alps::Parameters& parms, std::vector<std::pair<std::size_t, std::string> > const& qn);

template <class SymmGroup>
std::map<typename SymmGroup::charge,std::size_t> init_qn_charges
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states);

template <class SymmGroup>
std::map<alps::site_state<short>, std::pair<typename SymmGroup::charge, std::size_t> > init_coords
(std::vector<std::pair<std::size_t, std::string> > const & conserved_qn, alps::site_basis<short> const & states);


template <class SymmGroup>
class basis_converter {
public:
    typedef short I;
    typedef typename SymmGroup::charge charge_t;
    typedef std::size_t size_t;
    typedef std::pair<charge_t, size_t> coord_t;
            
    basis_converter (std::vector<std::pair<std::size_t, std::string> > const & conserved_qn,
                     alps::SiteBasisDescriptor<I> const & b)
    : states(alps::site_basis<I>(b))
    , charges(init_qn_charges<SymmGroup>(conserved_qn, states))
    , coords_(init_coords<SymmGroup>(conserved_qn, states))
    { }
    
    
    coord_t coords (alps::site_state<I> const & state) const
    {
        if (coords_.find(state) != coords_.end())
            return coords_.find(state)->second;
        else
            throw std::runtime_error("state not found!");
    }
    coord_t coords (size_t i) const { return coords(states[i]); }

    charge_t charge (alps::site_state<I> const & state) const
    {
        if (coords_.find(state) != coords_.end())
            return (coords_.find(state)->second).first;
        else
            throw std::runtime_error("state not found!");
    }
    charge_t charge (size_t i) const { return charge(states[i]); }

    size_t block_size (alps::site_state<I> const & state) const
    {
        charge_t c = charge(state);
        if (charges.find(c) != charges.end())
            return charges.find(c)->second;
        else
            throw std::runtime_error("state not found!");
    }
    size_t block_size (size_t i) const { return block_size(states[i]); }


    Index<SymmGroup> phys_dim () const
    {
        Index<SymmGroup> phys_i;
        
        for (typename std::map<charge_t,size_t>::const_iterator it = charges.begin();
             it != charges.end(); ++it)
            phys_i.insert( std::make_pair(it->first, it->second) );
        
        return phys_i;
    }
    
private:
    alps::site_basis<I> states;
    std::map<charge_t,size_t> charges;
    std::map<alps::site_state<I>, coord_t> coords_;
};

template <class Matrix, class SymmGroup>
class ALPSModel : public Model<Matrix, SymmGroup>
{
    typedef alps::SiteOperator SiteOperator;
    typedef alps::BondOperator BondOperator;
    
    typedef typename Matrix::value_type value_type;
    typedef boost::multi_array<value_type,2> alps_matrix;
    
    typedef typename basis_converter<SymmGroup>::I I;
    typedef alps::graph_helper<> graph_type;
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
public:
    typedef Hamiltonian<Matrix, SymmGroup> ham;        
    typedef typename ham::hamterm_t hamterm_t;        
    typedef typename ham::op_t op_t;
    
    typedef typename SymmGroup::charge charge;
    
    ALPSModel (const graph_type& lattice_, const alps::Parameters& parms_)
    : parms(parms_)
    , lattice(lattice_)
    , model(lattice, parms)
    {
        
        // Parsing conserved quantum numbers
        std::vector<std::string> tmp_qn;
        if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
            boost::char_separator<char> sep(" ,");
            std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
            tokenizer qn_tokens(qn_string, sep);
            for (tokenizer::iterator it=qn_tokens.begin();
                 it != qn_tokens.end();
                 it++)
            {
                tmp_qn.push_back(*it);
            }
        }
        
        // Load all possible basis
        for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
            alps::SiteBasisDescriptor<I> b = model.site_basis(type);
            alps::site_basis<I> states(b);
            // loop over states
            
            // TODO: QN only from type=0 vertex, what if there are more?
            if (type == 0)
                for (std::size_t i=0; i<b.size(); ++i)
                    if (std::find(tmp_qn.begin(), tmp_qn.end(), b[i].name()) != tmp_qn.end())
                        conserved_qn.push_back( std::make_pair(i, b[i].name()) );
            
            converter.push_back(basis_converter<SymmGroup>(conserved_qn, b));
            
            for (int i=0; i<states.size(); ++i) {
                charge c = converter[type].charge(i);
                size_t bsize = converter[type].block_size(i);
                // maquis::cout << "Inserting " << c << " for " << states[i] << std::endl;
                if (!tident[type].has_block(c, c))
                    tident[type].insert_block(Matrix::identity_matrix(bsize), c, c);
                int sign = (alps::is_fermionic(b, states[i])) ? -1 : 1;
                if (!tfill[type].has_block(c, c))
                    tfill[type].insert_block(Matrix::identity_matrix(bsize), c, c);
                tfill[type](converter[type].coords(i), converter[type].coords(i)) = sign;
            }
        }
        
        
        /*
         {
         maquis::cout << "BASIS:" << std::endl;
         alps::SiteBasisDescriptor<I> b = model.site_basis(0);
         alps::site_basis<I> states(b);
         for (typename std::map<int, std::vector<typename SymmGroup::charge> >::iterator it=tphys.begin();
         it != tphys.end();
         it++) {
         
         maquis::cout << "type " << it->first << ":" << std::endl;
         alps::SiteBasisDescriptor<I> b = model.site_basis(it->first);
         alps::site_basis<I> states(b);
         for (int i=0; i<it->second.size(); ++i) {
         maquis::cout << " " << i << ":" <<  " " << it->second[i] << " " << states[i] << std::endl;
         }
         
         }
         }
         */
        
        
        // site_term loop with cache to avoid recomputing matrices
        for (graph_type::site_iterator it=lattice.sites().first; it!=lattice.sites().second; ++it) {
            int p = lattice.vertex_index(*it);
            int type = lattice.site_type(*it);
            
            
            if (site_terms.find(type) == site_terms.end()) {
                typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator> > V;
                V  ops = model.site_term(type).template templated_split<value_type>();
                                                        
                for (int n=0; n<ops.size(); ++n) {
                    if (ops[n].get<0>().value() != 0.) {
                        SiteOperator op = ops[n].get<1>();
                        alps_matrix m = alps::get_matrix(value_type(), op, model.site_basis(type), parms, true);
                        
                        value_type coeff = ops[n].get<0>().value();
                        site_terms[type].push_back( coeff*convert_matrix(m, type) );
                    }
                    
                }
                
            }

            // Many site terms
            /*
             for (int n=0; n<site_terms[type].size(); ++n) {
             hamterm_t term;
             term.fill_operator = tident[type];
             term.operators.push_back( std::make_pair(p, site_terms[type][n]) );
             terms.push_back(term);
             }
             */
            
            // All site terms summed into one
            if (site_terms[type].size() > 0) {
                op_t op_matrix;
                for (int n=0; n<site_terms[type].size(); ++n)
                    op_matrix += site_terms[type][n];
                hamterm_t term;
                term.with_sign = false;
                term.fill_operator = tident[type];
                term.operators.push_back( std::make_pair(p, op_matrix) );
                terms.push_back(term);
            }
            
        }
        
        // bond_term loop
        for (graph_type::bond_iterator it=lattice.bonds().first; it!=lattice.bonds().second; ++it) {
            int p_s = lattice.source(*it);
            int p_t = lattice.target(*it);
            int type = lattice.bond_type(*it);
            int type_s = lattice.site_type(lattice.source(*it));
            int type_t = lattice.site_type(lattice.target(*it));
            
            bool wrap_pbc = boost::get(alps::boundary_crossing_t(),
                                       lattice.graph(),
                                       *it);
            
            BondOperator bondop = model.bond_term(type);
            
            typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
            alps::SiteBasisDescriptor<I> b1 = model.site_basis(type_s);
            alps::SiteBasisDescriptor<I> b2 = model.site_basis(type_t);
            
            
            V  ops = bondop.template templated_split<value_type>(b1,b2);
            for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                SiteOperator op1 = tit->get<1>();
                SiteOperator op2 = tit->get<2>();
                
                bool with_sign = fermionic(b1, op1, b2, op2);
                
                hamterm_t term;
                term.with_sign = with_sign;
                if (with_sign)
                    term.fill_operator = tfill[type_s];
                else
                    term.fill_operator = tident[type_s];
                {
                    alps_matrix m = alps::get_matrix(value_type(), op1, b1, parms, true);
                    value_type coeff = tit->get<0>().value();
                    op_t tmp;
                    if (with_sign && !wrap_pbc) {                            
                        gemm(tfill[type_s], convert_matrix(m, type_s), tmp); // Note inverse notation because of notation in operator.
                    } else
                        tmp = convert_matrix(m, type_s);
                    if (with_sign && wrap_pbc)
                        coeff *= -1.;
                    term.operators.push_back( std::make_pair(p_s, coeff*tmp) );
                }
                {
                    alps_matrix m = alps::get_matrix(value_type(), op2, b2, parms, true);
                    op_t tmp;
                    if (with_sign && wrap_pbc)
                        gemm(tfill[type_t], convert_matrix(m, type_t), tmp); // Note inverse notation because of notation in operator.
                    else
                        tmp = convert_matrix(m, type_t);
                    term.operators.push_back( std::make_pair(p_t, tmp) );
                }
                
                terms.push_back(term);
            }
            
        }
    }
            
    op_t const & get_identity() const
    {
        return tident[0];
    }
    
    Index<SymmGroup> get_phys() const
    {
        return converter[0].phys_dim();
    }
    
    Hamiltonian<Matrix, SymmGroup> H () const
    {
        return ham(get_phys(), get_identity(), terms);
    }
    
    typename SymmGroup::charge initc (BaseParameters& parms_) const
    {
        return init_charge<SymmGroup>(parms_, conserved_qn);
        /*        	typename SymmGroup::charge c = SymmGroup::IdentityCharge;
         for (int i=0; i<conserved_qn.size(); ++i) {
         if (conserved_qn.size() == 1)
         c = alps::evaluate<double>(static_cast<std::string>(parms[conserved_qn[0].second+"_total"]),parms)*2;
         else
         c[i] = alps::evaluate<double>(static_cast<std::string>(parms[conserved_qn[i].second+"_total"]),parms)*2;
         }
         return c;*/
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
                    charge c_i = converter[type].charge(i);
                    size_t bsize_i = converter[type].block_size(i);
                    charge c_j = converter[type].charge(j);
                    size_t bsize_j = converter[type].block_size(j);
                    
                    if (!newm.has_block(c_i, c_j))
                        newm.insert_block(Matrix(bsize_i, bsize_j, 0), c_i, c_j);
                    // Notation: going from state i to state j
                    newm(converter[type].coords(i), converter[type].coords(j)) = m[i][j];
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
    
    alps::Parameters const & parms;
    const graph_type& lattice;
    alps::model_helper<I> model;
    
    std::vector<basis_converter<SymmGroup> > converter;
    
    std::map<int, std::vector<op_t> > site_terms;
    mutable std::map<int, op_t> tident;
    mutable std::map<int, op_t> tfill;
    mutable std::map<int, std::vector<typename SymmGroup::charge> > tphys;
    std::vector<hamterm_t> terms;
    std::vector<std::pair<std::size_t, std::string> > conserved_qn;
    
};

namespace details {
    template <class T>
    int to_integer (alps::half_integer<T> const & qn_value)
    {
        return (qn_value.get_twice() % 2 == 0) ? alps::to_integer(qn_value) : qn_value.get_twice();
    }
}

// Loading Measurements
template <class Matrix, class SymmGroup>
Measurements<Matrix, SymmGroup> ALPSModel<Matrix, SymmGroup>::measurements () const
{
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    Measurements<Matrix, SymmGroup> meas;
    
    // TODO: hard coding site type = 0 !!
    int type = 0;
    
    meas.set_identity(tident[type]);
    
    {
        boost::regex expression("^MEASURE_AVERAGE\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            if (boost::regex_match(lhs, what, expression)) {
				alps::SiteBasisDescriptor<I> b = model.site_basis(type);

                if (model.has_bond_operator(it->value())) {
                	BondOperator bondop = model.get_bond_operator(it->value());

                    typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
                    V  ops = bondop.template templated_split<value_type>(b,b);
                    for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                        SiteOperator op1 = tit->get<1>();
                        SiteOperator op2 = tit->get<2>();

                        bool with_sign = fermionic(b, op1, b, op2);

                        mterm_t term;
                        term.type = mterm_t::Average;
                        std::ostringstream ss;
                        ss << what.str(1);
                        if (ops.size() > 1) ss << " (" << int(tit-ops.begin())+1 << ")";
                        term.name = ss.str();

                        term.fill_operator = (with_sign) ? tfill[type] : tident[type];
                        {
                        	alps_matrix m = alps::get_matrix(value_type(), op1, b, parms, true);
                        	op_t tmp;
                        	if (with_sign)
                        		gemm(tfill[type], convert_matrix(m, type), tmp); // Note inverse notation because of notation in operator.
                            else
                                tmp = convert_matrix(m, type);
                        	term.operators.push_back( std::make_pair(tit->get<0>().value()*tmp, b.is_fermionic(simplify_name(op1))) );
                        }
                        {
                            alps_matrix m = alps::get_matrix(value_type(), op2, b, parms, true);
                            term.operators.push_back( std::make_pair(convert_matrix(m, type), b.is_fermionic(simplify_name(op2))) );
                        }
	                    meas.add_term(term);
                    }
                } else {
                    mterm_t term;
                    term.type = mterm_t::Average;
                    term.name = what.str(1);

					SiteOperator op = make_site_term(it->value(), parms);
#ifndef NDEBUG
					if (b.is_fermionic(simplify_name(op)))
						throw std::runtime_error("Cannot measure local fermionic operators.");
#endif
					alps_matrix m = alps::get_matrix(value_type(), op, b, parms, true);

					term.operators.push_back( std::make_pair(convert_matrix(m, type), false) );
                    meas.add_term(term);
                }
            }
        }
    }
    
    {
        boost::regex expression("^MEASURE_LOCAL\\[(.*)]$");
        boost::smatch what;
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();
            if (boost::regex_match(lhs, what, expression)) {
				alps::SiteBasisDescriptor<I> b = model.site_basis(type);

                if (model.has_bond_operator(it->value())) {
                	BondOperator bondop = model.get_bond_operator(it->value());

                    typedef std::vector<boost::tuple<alps::expression::Term<value_type>,alps::SiteOperator,alps::SiteOperator > > V;
                    V  ops = bondop.template templated_split<value_type>(b,b);
                    for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                        SiteOperator op1 = tit->get<1>();
                        SiteOperator op2 = tit->get<2>();

                        bool with_sign = fermionic(b, op1, b, op2);

                        mterm_t term;
                        term.type = mterm_t::Local;
                        std::ostringstream ss;
                        ss << what.str(1);
                        if (ops.size() > 1) ss << " (" << int(tit-ops.begin())+1 << ")";
                        term.name = ss.str();

                        term.fill_operator = (with_sign) ? tfill[type] : tident[type];
                        {
                        	alps_matrix m = alps::get_matrix(value_type(), op1, b, parms, true);
                        	op_t tmp;
                        	if (with_sign)
                        		gemm(tfill[type], convert_matrix(m, type), tmp); // Note inverse notation because of notation in operator.
                            else
                                tmp = convert_matrix(m, type);
                        	term.operators.push_back( std::make_pair(tit->get<0>().value()*tmp, b.is_fermionic(simplify_name(op1))) );
                        }
                        {
                            alps_matrix m = alps::get_matrix(value_type(), op2, b, parms, true);
                            term.operators.push_back( std::make_pair(convert_matrix(m, type), b.is_fermionic(simplify_name(op2))) );
                        }
	                    meas.add_term(term);
                    }
                } else {
                    mterm_t term;
                    term.type = mterm_t::Local;
                    term.name = what.str(1);

					SiteOperator op = make_site_term(it->value(), parms);
#ifndef NDEBUG
					if (b.is_fermionic(simplify_name(op)))
						throw std::runtime_error("Cannot measure local fermionic operators.");
#endif

					alps_matrix m = alps::get_matrix(value_type(), op, b, parms, true);

					term.operators.push_back( std::make_pair(convert_matrix(m, type), false) );
                    meas.add_term(term);
                }
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
            term.fill_operator = tident[type];
            
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
                    alps_matrix m = alps::get_matrix(value_type(), op, b, parms, true);
                    bool f = b.is_fermionic(simplify_name(op));
                    term.operators.push_back( std::make_pair(convert_matrix(m, type), f) );
                    if (f) ++f_ops;
                }
                if (term.operators.size() == 1) {
                    term.operators.push_back(term.operators[0]);
                    if (term.operators[1].second) ++f_ops;
                }
                
                if (f_ops > 0) {
                    term.fill_operator = tfill[type];
                }
                
#ifndef NDEBUG
                if (f_ops % 2 != 0)
                    throw std::runtime_error("Number of fermionic operators has to be even.");
#endif
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
            term.fill_operator = tident[type];
            
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
                
                boost::char_separator<char> sep(":");
                tokenizer corr_tokens(value, sep);
                for (tokenizer::iterator it2=corr_tokens.begin();
                     it2 != corr_tokens.end();
                     it2++)
                {
                    SiteOperator op = make_site_term(*it2, parms);
                    alps_matrix m = alps::get_matrix(value_type(), op, b, parms, true);
                    bool f = b.is_fermionic(simplify_name(op));
                    term.operators.push_back( std::make_pair(convert_matrix(m, type), f) );
                    if (f) ++f_ops;
                }
                if (term.operators.size() == 1) {
                    term.operators.push_back(term.operators[0]);
                    if (term.operators[1].second) ++f_ops;
                }
                
                
                if (f_ops > 0) {
                    term.fill_operator = tfill[type];
                }

#ifndef NDEBUG
                if (f_ops % 2 != 0)
                    throw std::runtime_error("Number of fermionic operators has to be even.");
#endif
                meas.add_term(term);
            }
        }
    }
    
    
    return meas;
}

#endif
