#ifndef APP_ALPS_MODEL_H
#define APP_ALPS_MODEL_H

#include "hamiltonian.h"

#include <alps/parameter.h>
#include <alps/lattice.h>
#include <alps/model.h>

namespace app {
    
    template <class SymmGroup>
    typename SymmGroup::charge convert_alps (alps::site_state<short> const & state);

    template <>
    U1::charge convert_alps<U1> (alps::site_state<short> const & state)
    {
        assert(state.size() == 1);
        return get_quantumnumber(state, 0).get_twice();
    }
    template <>
    TwoU1::charge convert_alps<TwoU1> (alps::site_state<short> const & state)
    {
        assert(state.size() == 2);
        TwoU1::charge ret;
        ret[0] = get_quantumnumber(state, 0).get_twice();
        ret[1] = get_quantumnumber(state, 1).get_twice();
        return ret;
    }

    
    template <class Matrix, class SymmGroup>
    class ALPSModel : public Hamiltonian<Matrix, SymmGroup>
    {
        typedef alps::SiteOperator SiteOperator;
        typedef alps::BondOperator BondOperator;
        typedef boost::multi_array<double,2> alps_matrix;
        typedef short I;
        typedef alps::graph_helper<> graph_type;
        
    public:
        typedef Hamiltonian_Term<Matrix, SymmGroup> hamterm_t;
        typedef typename hamterm_t::op_t op_t;
        
        typedef typename SymmGroup::charge charge;
        
        ALPSModel (const graph_type& lattice, const alps::Parameters& parms)
        {
            
            alps::model_helper<I> model(lattice, parms);
            
            // Load all possible basis
            for (int type=0; type<=alps::maximum_vertex_type(lattice.graph()); ++type) {
                alps::site_basis<I> states = model.site_basis(type);
                // loop over states
                for (int i=0; i<states.size(); ++i) {
                    charge c = convert_alps<SymmGroup>(states[i]);
                    tphys[type].insert(std::make_pair(c, 1));
                    tident[type].insert_block(Matrix(1, 1, 1), c, c);
                }
            }

            /*
            std::cout << "BASIS:" << std::endl;
            for (typename std::map<int, Index<SymmGroup> >::iterator it=tphys.begin();
                 it != tphys.end();
                 it++) {
                
                std::cout << it->first << " " << it->second << std::endl;
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
                    V  ops = model.site_term(type).templated_split<double>();
#ifndef NDEBUG
                    if (ops.size() != 1) {
                        std::runtime_error("SiteOperator not of length one.");
                    }
#endif
                    if (ops[0].get<0>().value() != 0.) {
                        SiteOperator op = ops[0].get<1>();
                        alps_matrix m = alps::get_matrix(double(),
                                                         op,
                                                         model.site_basis(type),
                                                         parms,
                                                         true);
                        
                        for (int i=0; i<m.shape()[0]; ++i) {
                            for (int j=0; j<m.shape()[1]; ++j) {
                                newm.insert_block(Matrix(1, 1, ops[0].get<0>().value()*m[i][j]),
                                                  tphys[type][i].first,
                                                  tphys[type][j].first);
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
                    term.fill_operator = tident[0];
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
                
                BondOperator bondop = model.bond_term(type);
                
                typedef std::vector<boost::tuple<alps::expression::Term<double>,alps::SiteOperator,alps::SiteOperator > > V;
                alps::SiteBasisDescriptor<I> b1 = model.site_basis(type_s);
                alps::SiteBasisDescriptor<I> b2 = model.site_basis(type_t);
                
                hamterm_t term;
                term.fill_operator = tident[0];
                
                V  ops = bondop.template templated_split<double>(b1,b2);
                for (typename V::iterator tit=ops.begin(); tit!=ops.end();++tit) {
                    SiteOperator op1 = tit->get<1>();
                    SiteOperator op2 = tit->get<2>();
                    
                    {
                        alps_matrix m = alps::get_matrix(double(), op1, b1, parms, true);
                        op_t newm;
                        for (int i=0; i<m.shape()[0]; ++i) {
                            for (int j=0; j<m.shape()[1]; ++j) {
                                newm.insert_block(Matrix(1, 1, m[i][j]), tphys[type_s][i].first, tphys[type_s][j].first);
                            }
                        }
                        term.operators.push_back( std::make_pair(p_s, tit->get<0>().value()*newm) );
                    }
                    {
                        alps_matrix m = alps::get_matrix(double(), op2, b2, parms, true);
                        op_t newm;
                        for (int i=0; i<m.shape()[0]; ++i) {
                            for (int j=0; j<m.shape()[1]; ++j) {
                                newm.insert_block(Matrix(1, 1, m[i][j]), tphys[type_t][i].first, tphys[type_t][j].first);
                            }
                        }
                        term.operators.push_back( std::make_pair(p_t, newm) );
                    }
                    
                }

                terms.push_back(term);
                
            }
        }
        
        int n_terms(TermsType what) const
        {
            return terms.size();
        }
        hamterm_t operator[](int i) const
        {
            return terms[i];
        }
        
        op_t get_identity() const
        {
            return tident[0];
        }
        
        Index<SymmGroup> get_phys() const
        {
            // TODO: do we need multiple site types?
            return tphys[0];
        }
        
    private:
        
        charge create_charge (alps::site_state<I> const & state) const;
        
        
        std::vector<hamterm_t> terms;
        std::map<int, std::pair<bool, op_t> > site_terms;
        mutable std::map<int, op_t> tident;
        mutable std::map<int, Index<SymmGroup> > tphys;
        
    };
    
    
} // namespace

#endif