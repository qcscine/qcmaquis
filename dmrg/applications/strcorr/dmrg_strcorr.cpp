
#include <iterator>
#include <iostream>
#include <sys/stat.h>

#include "dmrg/utils/DmrgParameters2.h"
#include "utils/timings.h"

#include "dmrg_strcorr.h"


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
    
    if (model["MODEL"] != "optical_lattice")
        throw std::runtime_error("This application works only with `optical_lattice` continuum models.");
    
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);

    try {
        
        StrCorr sim(parms, model);
        
        for (int l=1; l<=8; ++l) {
            maquis::cout << "Measure single-site string operator, size " << l*model["Ndiscr"] << "." << std::endl;
            sim.measure_ss_string_unit(l*model["Ndiscr"]);
            maquis::cout << "Measure unit-cell string operator, size " << l << "." << std::endl;
            sim.measure_uc_string(l*model["Ndiscr"]);
        }
        
    } catch (std::exception & e) {
        maquis::cerr << "Exception thrown!" << std::endl;
        maquis::cerr << e.what() << std::endl;
        exit(1);
    }
    
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
}

