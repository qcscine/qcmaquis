#include "dmrg_tevol.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "utils/data_collector.hpp"
#include "dmrg/utils/DmrgParameters2.h"
#include "utils/timings.h"

#include <boost/function.hpp>
#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>

std::string guess_alps_symmetry(ModelParameters& model)
{
    std::map<int, std::string> symm_names;
    symm_names[0] = "none";
    symm_names[1] = "u1";
    symm_names[2] = "2u1";
    
    int n=0;
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    if (model.defined("CONSERVED_QUANTUMNUMBERS")) {
        boost::char_separator<char> sep(" ,");
        std::string qn_string = model["CONSERVED_QUANTUMNUMBERS"];
        tokenizer qn_tokens(qn_string, sep);
        for (tokenizer::iterator it=qn_tokens.begin(); it != qn_tokens.end(); it++) {
            if (model.defined(*it + "_total"))
                n += 1;
        }
    }
    return symm_names[n];
}

// factory
void dmrg_tevol(DmrgParameters & parms, ModelParameters & model)
{
    std::map<std::string, boost::function<void (DmrgParameters & p, ModelParameters & m)> > factory_map;
    
    bool img_only = ( parms["nsweeps"] == parms["nsweeps_img"] );
    
#ifdef HAVE_TrivialGroup
    if (img_only)
        factory_map["none"] = run_tevol<matrix, TrivialGroup>;
    else
        factory_map["none"] = run_tevol<cmatrix, TrivialGroup>;
#endif
#ifdef HAVE_U1
    if (img_only)
        factory_map["u1"] = run_tevol<matrix, U1>;
    else
        factory_map["u1"] = run_tevol<cmatrix, U1>;
#endif
#ifdef HAVE_TwoU1
    if (img_only)
        factory_map["2u1"] = run_tevol<matrix, TwoU1>;
    else
        factory_map["2u1"] = run_tevol<cmatrix, TwoU1>;
#endif
    
    std::string symm_name;
    if (parms["model_library"] == "alps")
        symm_name = guess_alps_symmetry(model);
    else
        symm_name = parms["symmetry"].str();
    
    if (factory_map.find(symm_name) != factory_map.end())
        factory_map[symm_name](parms, model);
    else
        throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");    
}

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
        exit(1);
    }
    
    maquis::cout.precision(10);
    
    std::ifstream param_file(argv[1]);
    if (!param_file) {
        maquis::cerr << "Could not open parameter file." << std::endl;
        exit(1);
    }
    DmrgParameters parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        maquis::cerr << "Could not open model file." << std::endl;
        exit(1);
    }
    ModelParameters model(model_file);
    
    DCOLLECTOR_SET_SIZE(gemm_collector, parms["max_bond_dimension"]+1)
    DCOLLECTOR_SET_SIZE(svd_collector, parms["max_bond_dimension"]+1)
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    try {
        dmrg_tevol(parms, model);
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
    DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
    
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
}

