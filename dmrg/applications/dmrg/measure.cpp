#include "measure.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

#include "dmrg/utils/DmrgParameters2.h"


// factory
void dmrg(DmrgParameters & parms, ModelParameters & model)
{
    std::map<std::string, boost::function<void (DmrgParameters & p, ModelParameters & m)> > factory_map;
    
    maquis::cout << "This binary contains symmetries: ";
#ifdef HAVE_TrivialGroup
    factory_map["none"] = run_dmrg<TrivialGroup>;
    maquis::cout << "none ";
#endif
#ifdef HAVE_U1
    factory_map["u1"] = run_dmrg<U1>;
    maquis::cout << "u1 ";
#endif
#ifdef HAVE_TwoU1
    factory_map["2u1"] = run_dmrg<TwoU1>;
    maquis::cout << "2u1 ";
#endif
    maquis::cout << std::endl;

    
    std::string symm_name = parms.get<std::string>("symmetry");
    
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
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    
    try {
        dmrg(parms, model);
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }    
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
}

