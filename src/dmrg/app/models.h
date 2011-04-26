#ifndef APP_MODELS_H
#define APP_MODELS_H

#include "hamiltonian.h"
#include "lattice.h"

namespace app {

template <class Matrix, class SymmGroup>
void model_parser (std::string const & type, std::string const & fname,
                   Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H, typename SymmGroup::charge& initc);
}


/*
#include "utils/DmrgParameters.h"
#include <alps/parameter.h>
#include "alps_lattice.hpp"
#include "hamiltonians.hpp"
#include "alps_model.hpp"

#include <boost/tokenizer.hpp>

    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    
    template<class Matrix, class SymmGroup>
    Hamiltonian<Matrix, SymmGroup> * hamil_factory(BaseParameters & model, Lattice const & lattice)
    {
#ifdef UseTwoU1
        if (model.get<std::string>("model") == std::string("fermi_hubbard"))
            return new TwoU1_FermiHubbard<Matrix>(lattice, model);
        else {
            throw std::runtime_error("Don't know this model!");
            return NULL;
        }
#else
        if (model.get<std::string>("model") == std::string("heisenberg"))
            return new Heisenberg<Matrix>(lattice, model.get<double>("Jxy"), model.get<double>("Jz"));
        else if (model.get<std::string>("model") == std::string("HCB"))
            return new HCB<Matrix>(lattice);
        else if (model.get<std::string>("model") == std::string("FreeFermions"))
            return new FreeFermions<Matrix>(lattice, model.get<double>("t"));
        else {
            throw std::runtime_error("Don't know this model!");
            return NULL;
        }
#endif
    }


    template <class SymmGroup>
    typename SymmGroup::charge init_qn (alps::Parameters const& parms);
    
    template <class Matrix, class SymmGroup>
    void model_parser (std::string const & type, std::string const & fname,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H, typename SymmGroup::charge& initc)
    {
        
        if (type == "alps") {
            std::ifstream ifs(fname.c_str());
            alps::Parameters parms(ifs);
            ALPSLattice * lat_tmp = new ALPSLattice(parms);
            ALPSModel<Matrix, SymmGroup> * Htmp = new ALPSModel<Matrix, SymmGroup>(lat_tmp->alps_graph(), parms);
            lattice = lat_tmp;
            H = Htmp;
            initc = Htmp->init_qn(parms);
        } else if (type == "coded") {
            std::ifstream model_file(fname.c_str());
            std::ifstream ifs(fname.c_str());
            alps::Parameters parms(ifs);
            lattice = new ALPSLattice(parms);
            ModelParameters model(model_file);
            H = hamil_factory<Matrix, SymmGroup>(model, *lattice);
#ifdef UseTwoU1
            initc[0] = model.get<int>("u1_total_charge1");
            initc[1] = model.get<int>("u1_total_charge2");
#else
            initc = model.get<int>("u1_total_charge");
#endif
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
*/

#endif
