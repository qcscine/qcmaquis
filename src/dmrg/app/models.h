#ifndef APP_MODELS_H
#define APP_MODELS_H

#include <fstream>

#include "hamiltonian.h"
#include "lattice.h"

#include <alps/parameter.h>
#include "alps_lattice.hpp"
#include "hamiltonians.hpp"
#include "alps_model.hpp"

#include <boost/tokenizer.hpp>

namespace app {
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    
    //    void model_parser (std::string const & type, std::ifstream & ifs,
    //                        Lattice & lattice, Hamiltonian<> & H);
    
    template <class SymmGroup>
    typename SymmGroup::charge init_qn (alps::Parameters const& parms);
    
    template <class Matrix, class SymmGroup>
    void model_parser (std::string const & type, std::ifstream & ifs,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H, typename SymmGroup::charge& initc)
    {
        
        if (type == "alps") {
            alps::Parameters parms(ifs);
            ALPSLattice * tmp = new ALPSLattice(parms);
            H = new ALPSModel<Matrix, SymmGroup>(tmp->alps_graph(), parms);
            lattice = tmp;
            initc = init_qn<SymmGroup>(parms);
        } else {
            throw std::runtime_error("Don't know this type of lattice / model!");
        }
        
    }
    
    template <>
    U1::charge init_qn<U1> (alps::Parameters const& parms)
    {
        U1::charge ret = U1::SingletCharge;
        
        boost::char_separator<char> sep(" ,");
        
        if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
            std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
            tokenizer qn_tokens(qn_string, sep);
            int count = 0;
            for (tokenizer::iterator it=qn_tokens.begin();
                 it != qn_tokens.end();
                 it++)
            {
                if (parms.defined(*it+"_total")) {
                    assert( count<1 );
                    ret = alps::evaluate<double>(static_cast<std::string>(parms[*it+"_total"]),parms)*2;
                    count++;
                }
            }
        }
        
        return ret;
    }

    template <>
    TwoU1::charge init_qn<TwoU1> (alps::Parameters const& parms)
    {
        TwoU1::charge ret = TwoU1::SingletCharge;
        
        boost::char_separator<char> sep(" ,");
        
        if (parms.defined("CONSERVED_QUANTUMNUMBERS")) {
            std::string qn_string = parms["CONSERVED_QUANTUMNUMBERS"];
            tokenizer qn_tokens(qn_string, sep);
            int count = 0;
            for (tokenizer::iterator it=qn_tokens.begin();
                 it != qn_tokens.end();
                 it++)
            {
                if (parms.defined(*it+"_total")) {
                    assert( count<2 );
                    ret[count] = alps::evaluate<double>(static_cast<std::string>(parms[*it+"_total"]),parms)*2;
                    count++;
                }
            }
        }
        
        return ret;
    }
    
}

#endif
