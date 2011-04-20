#ifndef B_DMRGPARAMETERS_H
#define B_DMRGPARAMETERS_H

#include <boost/program_options.hpp>
#include <boost/tokenizer.hpp>
#include <boost/lexical_cast.hpp>

#include <boost/any.hpp>

#include <fstream>
#include <iostream>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

#include "utils/DmrgParameters.h"

class b_DmrgParameters : public BaseParameters
{
public:
    b_DmrgParameters(std::ifstream& param_file)
    {
        using namespace boost::program_options;
        
        try {
            add_option("truncation_initial",value<double>(),"Initial value for the truncation error");
            add_option("truncation_final",value<double>(),"Final value for the truncation");
            
            add_option("init_bond_dimension",value<std::size_t>()->default_value(5),"");
            add_option("max_bond_dimension",value<std::size_t>(),"");
            
            add_option("alpha_initial",value<double>()->default_value(1e-2),"");
            add_option("alpha_main",value<double>()->default_value(1e-6),"");
            add_option("alpha_final",value<double>()->default_value(1e-20),"");
            
            add_option("eigensolver",value<std::string>()->default_value(std::string("ARPACK")),"");
            add_option("arpack_tol",value<double>()->default_value(1e-8),"");
            add_option("arpack_ncv",value<int>()->default_value(20),"");
            add_option("ietl_jcd_tol",value<double>()->default_value(1e-8),"");
            add_option("ietl_jcd_gmres",value<int>()->default_value(5),"");
            add_option("ietl_jcd_maxiter",value<int>()->default_value(100),"");
            
            add_option("nsweeps",value<int>(),"");
            add_option("nmainsweeps",value<int>(),"");
            add_option("ngrowsweeps",value<int>(),"");
            
            add_option("resultfile",value<std::string>(),"");
            add_option("chkpfile",value<std::string>(),"");
            add_option("initfile",value<std::string>()->default_value(std::string()),"");
            
            add_option("donotsave",value<int>()->default_value(0),"");
            add_option("run_seconds",value<int>()->default_value(0),"");
            add_option("storagedir",value<std::string>()->default_value(std::string()),"");
            add_option("use_compressed",value<int>()->default_value(0),"");
            add_option("calc_h2",value<int>()->default_value(0),"");
            add_option("seed",value<int>()->default_value(42),"");
            
            add_option("init_state", value<std::string>()->default_value("default"),"");
            
            store(parse_config_file(param_file, config), vm);
            notify(vm);
        } catch(std::exception &e) {
            std::cerr << "Error reading parameter file." << std::endl;
            throw e;
        }
    }
};
    
class b_ModelParameters : public BaseParameters
{
public:
    b_ModelParameters(std::ifstream& param_file)
    {
        using namespace boost::program_options;
        
        try {
            add_option("model", value<std::string>(), "");
            add_option("lattice", value<std::string>(), "");
            add_option("alps_lattice",value<std::string>(),"");
            
            add_option("L", value<int>(), "");
            add_option("W", value<int>(), "");
            
            add_option("Jxy", value<double>(), "");
            add_option("Jz", value<double>(), "");
            
            add_option("U", value<double>(), "");
            add_option("t", value<double>()->default_value(1.), "");
            add_option("t1", value<double>()->default_value(1.), "");
            add_option("t2", value<double>()->default_value(1.), "");
            
            add_option("theta", value<double>(), "");
            add_option("h0", value<double>()->default_value(0), "");
            add_option("pin", value<int>()->default_value(3), "");
            
            add_option("K0", value<double>()->default_value(1), "");
            add_option("K1", value<double>()->default_value(0), "");
            
            add_option("penalty", value<double>()->default_value(10), "");
            add_option("twist", value<double>()->default_value(0), "");
            add_option("move", value<double>()->default_value(0), "");
            
            add_option("u1_total_charge", value<int>()->default_value(0), "");
            add_option("u1_total_charge1", value<int>()->default_value(0), "");
            add_option("u1_total_charge2", value<int>()->default_value(0), "");
            
            store(parse_config_file(param_file, config), vm);
            notify(vm);
        } catch(std::exception &e) {
            std::cerr << "Error reading parameter file." << std::endl;
            std::cerr << e.what() << endl;
            throw e;
        }
    }
};

#endif
