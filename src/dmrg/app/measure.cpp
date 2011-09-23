#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "utils/DmrgParameters2.h"
#include "measure.h"

// factory
void dmrg(DmrgParameters & parms, ModelParameters & model)
{
    if (parms.get<std::string>("symmetry") == "null")
        run_dmrg<NullGroup>(parms, model);
    else if (parms.get<std::string>("symmetry") == "u1")
        run_dmrg<U1>(parms, model);
    else if (parms.get<std::string>("symmetry") == "2u1")
        run_dmrg<TwoU1>(parms, model);
    else
        throw std::runtime_error("Don't know this symmetry group.");    
}

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        cout << "Usage: <parms> <model_parms>" << endl;
        exit(1);
    }
    
    cout.precision(10);

    std::ifstream param_file(argv[1]);
    if (!param_file) {
        cerr << "Could not open parameter file." << endl;
        exit(1);
    }
    DmrgParameters parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        cerr << "Could not open model file." << endl;
        exit(1);
    }
    ModelParameters model(model_file);
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    
    dmrg(parms, model);
    
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    cout << "Task took " << elapsed << " seconds." << endl;
}

