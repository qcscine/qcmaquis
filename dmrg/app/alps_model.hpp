#ifndef APP_ALPS_MODEL_H
#define APP_ALPS_MODEL_H

#include "hamiltonian.h"
#include "measurements.h"

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

namespace app {
    
    template <class SymmGroup>
    typename SymmGroup::charge convert_alps (alps::site_state<short> const & state, std::vector<std::pair<int, std::string> > const& qn);
    
    template <class SymmGroup>
    typename SymmGroup::charge init_charge (const alps::Parameters& parms, std::vector<std::pair<int, std::string> > const& qn);
    
    template <class Matrix, class SymmGroup>
    class ALPSModel
    {
        typedef alps::SiteOperator SiteOperator;
        typedef alps::BondOperator BondOperator;
        typedef boost::multi_array<double,2> alps_matrix;
        typedef short I;
        typedef alps::graph_helper<> graph_type;
        
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
        
    public:
        typedef Hamiltonian<Matrix, SymmGroup> ham;        
        typedef typename ham::hamterm_t hamterm_t;        
        typedef typename ham::op_t op_t;
        
        typedef typename SymmGroup::charge charge;
        
        ALPSModel (const graph_type& lattice_, const alps::Parameters& parms) :
		lattice(lattice_),
		model(lattice, parms)
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
#ifdef NDEBUG
            }
#else
            } else {
                std::runtime_error("No conserved quantum numbers defined!");
            }
#endif
            
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
                            
                            for (int i=0; i<states.size(); ++i) {
                                charge c = convert_alps<SymmGroup>(states[i], conserved_qn);
                                //std::cout << "Inserting " << c << " for " << states[i] << std::endl;
                                tphys[type].push_back(c);
                                tident[type].insert_block(Matrix(1, 1, 1), c, c);
                                int sign = (alps::is_fermionic(b, states[i])) ? -1 : 1;
                                tfill[type].insert_block(Matrix(1, 1, sign), c, c);
                            }
            }
            
            
            /*
             {
             std::cout << "BASIS:" << std::endl;
             alps::SiteBasisDescriptor<I> b = model.site_basis(0);
             alps::site_basis<I> states(b);
             for (typename std::map<int, std::vector<typename SymmGroup::charge> >::iterator it=tphys.begin();
             it != tphys.end();
             it++) {
             
             std::cout << "type " << it->first << ":" << std::endl;
             alps::SiteBasisDescriptor<I> b = model.site_basis(it->first);
             alps::site_basis<I> states(b);
             for (int i=0; i<it->second.size(); ++i) {
             std::cout << " " << i << ":" <<  " " << it->second[i] << " " << states[i] << std::endl;
             }
             
             }
             }
             */
            
            
            // site_term loop with cache to avoid recomputing matrices
            for (graph_type::site_iterator it=lattice.sites().first; it!=lattice.sites().second; ++it) {
                int p = lattice.vertex_index(*it);
                int type = lattice.site_type(*it);
                
                op_t newm;
                bool used;
                
                if (site_terms.find(type) == site_terms.end()) {
                    typedef std::vector<boost::tuple<alps::expression::Term<double>,alps::SiteOperator> > V;
                    V  ops = model.site_term(type).template templated_split<double>();
#ifndef NDEBUG
                    if (ops.size() != 1) {
                        std::runtime_error("SiteOperator not of length one.");
                    }
#endif
                    if (ops[0].get<0>().value() != 0.) {
                        SiteOperator op = ops[0].get<1>();
                        alps_matrix m = alps::get_matrix(double(), op, model.site_basis(type), parms, true);
                        
                        for (int i=0; i<m.shape()[0]; ++i) {
                            for (int j=0; j<m.shape()[1]; ++j) {
                                if (m[i][j] != 0.)
                                    // Notation: going from state i to state j
                                    newm.insert_block(Matrix(1, 1, ops[0].get<0>().value()*m[i][j]),
                                                      tphys[type][i],
                                                      tphys[type][j]);
                                    }
                        }
                        used = true;
                    } else {
                        used = false;
                    }
                    
                    site_terms[type] = std::make_pair(used, newm);
                } else {
                    newm = site_terms[type].second;
                    used = site_terms[type].first;
                }
                
                if (used) {
                    hamterm_t term;
                    term.fill_operator = tident[type];
                    term.operators.push_back( std::make_pair(p, newm) );
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
                
                typedef std::vector<boost::tuple<alps::expression::Term<double>,alps::SiteOperator,alps::SiteOperator > > V;
                alps::SiteBasisDescriptor<I> b1 = model.site_basis(type_s);
                alps::SiteBasisDescriptor<I> b2 = model.site_basis(type_t);
                
                
                V  ops = bondop.template templated_split<double>(b1,b2);
                for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                    SiteOperator op1 = tit->get<1>();
                    SiteOperator op2 = tit->get<2>();
                    
                    bool with_sign = fermionic(b1, op1, b2, op2);
                    
                    hamterm_t term;
                    if (with_sign)
                        term.fill_operator = tfill[type_s];
                    else
                        term.fill_operator = tident[type_s];
                    {
                        alps_matrix m = alps::get_matrix(double(), op1, b1, parms, true);
                        double coeff = tit->get<0>().value();
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
                        alps_matrix m = alps::get_matrix(double(), op2, b2, parms, true);
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
            Index<SymmGroup> phys;
            // TODO: do we need multiple site types?
            for (int i=0; i<tphys[0].size(); ++i) {
                phys.insert( std::make_pair(tphys[0][i], 1) );
            }
            return phys;
        }
        
        Hamiltonian<Matrix, SymmGroup> get_hamiltonian () const
        {
            return ham(get_phys(), get_identity(), terms);
        }
        
        typename SymmGroup::charge init_qn (const alps::Parameters& parms) const
        {
            return init_charge<SymmGroup>(parms, conserved_qn);
            /*        	typename SymmGroup::charge c = SymmGroup::SingletCharge;
             for (int i=0; i<conserved_qn.size(); ++i) {
             if (conserved_qn.size() == 1)
             c = alps::evaluate<double>(static_cast<std::string>(parms[conserved_qn[0].second+"_total"]),parms)*2;
             else
             c[i] = alps::evaluate<double>(static_cast<std::string>(parms[conserved_qn[i].second+"_total"]),parms)*2;
             }
             return c;*/
        }
        
        Measurements<Matrix, SymmGroup> parse_measurements (alps::Parameters const & parms) const;
        
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
                    if (m[i][j] != 0.)
                        // Notation: going from state i to state j
                        newm.insert_block(Matrix(1, 1, m[i][j]), tphys[type][i], tphys[type][j]);
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
        
        const graph_type& lattice;
        alps::model_helper<I> model;
        
        std::map<int, std::pair<bool, op_t> > site_terms;
        mutable std::map<int, op_t> tident;
        mutable std::map<int, op_t> tfill;
        mutable std::map<int, std::vector<typename SymmGroup::charge> > tphys;
        std::vector<hamterm_t> terms;
        std::vector<std::pair<int, std::string> > conserved_qn;
        
    };


    // Symmetry dependent implementation

    // U1 Symmetry
    template <>
    U1::charge init_charge<U1> (const alps::Parameters& parms, std::vector<std::pair<int, std::string> > const& qn)
    {
        assert(qn.size() == 1);
        U1::charge c = U1::SingletCharge;
        if (parms.defined(qn[0].second+"_total")) {
            c = alps::evaluate<double>(static_cast<std::string>(parms[qn[0].second+"_total"]),parms)*2;
        }
        return c;
    }

    template <>
    U1::charge convert_alps<U1> (alps::site_state<short> const & state, std::vector<std::pair<int, std::string> > const& qn)
    {
        assert(qn.size() == 1);
        return get_quantumnumber(state, qn[0].first).get_twice();
    }

    // TwoU1 Symmetry
    template <>
    TwoU1::charge convert_alps<TwoU1> (alps::site_state<short> const & state, std::vector<std::pair<int, std::string> > const& qn)
    {
        assert(qn.size() == 2);
        TwoU1::charge ret;
        ret[0] = get_quantumnumber(state, qn[0].first).get_twice();
        ret[1] = get_quantumnumber(state, qn[1].first).get_twice();
        return ret;
    }

    template <>
    TwoU1::charge init_charge<TwoU1> (const alps::Parameters& parms, std::vector<std::pair<int, std::string> > const& qn)
    {
        assert(qn.size() == 2);
        TwoU1::charge c = TwoU1::SingletCharge;
        if (parms.defined(qn[0].second+"_total")) {
            c[0] = alps::evaluate<double>(static_cast<std::string>(parms[qn[0].second+"_total"]),parms)*2;
        }
        if (parms.defined(qn[1].second+"_total")) {
            c[1] = alps::evaluate<double>(static_cast<std::string>(parms[qn[1].second+"_total"]),parms)*2;
        }
        return c;
    }


    // Loading Measurements
    template <class Matrix, class SymmGroup>
    Measurements<Matrix, SymmGroup> ALPSModel<Matrix, SymmGroup>::parse_measurements (alps::Parameters const & parms) const
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

                        typedef std::vector<boost::tuple<alps::expression::Term<double>,alps::SiteOperator,alps::SiteOperator > > V;
                        V  ops = bondop.template templated_split<double>(b,b);
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
                            	alps_matrix m = alps::get_matrix(double(), op1, b, parms, true);
                            	op_t tmp;
                            	if (with_sign)
                            		gemm(tfill[type], convert_matrix(m, type), tmp); // Note inverse notation because of notation in operator.
                                else
                                    tmp = convert_matrix(m, type);
                            	term.operators.push_back( std::make_pair(tit->get<0>().value()*tmp, b.is_fermionic(simplify_name(op1))) );
                            }
                            {
                                alps_matrix m = alps::get_matrix(double(), op2, b, parms, true);
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
							std::runtime_error("Cannot measure local fermionic operators.");
#endif
						alps_matrix m = alps::get_matrix(double(), op, b, parms, true);

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

                        typedef std::vector<boost::tuple<alps::expression::Term<double>,alps::SiteOperator,alps::SiteOperator > > V;
                        V  ops = bondop.template templated_split<double>(b,b);
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
                            	alps_matrix m = alps::get_matrix(double(), op1, b, parms, true);
                            	op_t tmp;
                            	if (with_sign)
                            		gemm(tfill[type], convert_matrix(m, type), tmp); // Note inverse notation because of notation in operator.
                                else
                                    tmp = convert_matrix(m, type);
                            	term.operators.push_back( std::make_pair(tit->get<0>().value()*tmp, b.is_fermionic(simplify_name(op1))) );
                            }
                            {
                                alps_matrix m = alps::get_matrix(double(), op2, b, parms, true);
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
							std::runtime_error("Cannot measure local fermionic operators.");
#endif
						alps_matrix m = alps::get_matrix(double(), op, b, parms, true);

						term.operators.push_back( std::make_pair(convert_matrix(m, type), false) );
	                    meas.add_term(term);
                    }
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
                term.fill_operator = tfill[type];
                
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
                        alps_matrix m = alps::get_matrix(double(), op, b, parms, true);
                        bool f = b.is_fermionic(simplify_name(op));
                        term.operators.push_back( std::make_pair(convert_matrix(m, type), f) );
                        if (f) ++f_ops;
                    }
                    if (term.operators.size() == 1) {
                        term.operators.push_back(term.operators[0]);
                        if (term.operators[1].second) ++f_ops;
                    }
                    
#ifndef NDEBUG
                    if (f_ops % 2 != 0)
                        std::runtime_error("Number of fermionic operators has to be even.");
#endif
                    meas.add_term(term);
                }
            }
        }
        
        
        return meas;
    }

} // namespace

#endif
