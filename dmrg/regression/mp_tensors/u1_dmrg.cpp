#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double > Matrix;


#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/mp_tensors/optimize.h"
#include "dmrg/deprecated/mpos/measurements.h"

#include "dmrg/deprecated/mpos/adjacency.h"
#include "dmrg/deprecated/mpos/generate_mpo.h"
#include "dmrg/deprecated/mpos/hamiltonians.h"

#include "dmrg/utils/DmrgParameters.h"
#include "utils/timings.h"

#include "dmrg/utils/storage.h"

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

adj::Adjacency * adj_factory(ModelParameters & model)
{
    if (model["lattice"] == std::string("square_lattice"))
        return new adj::SquareAdj(model["L"], model["W"]);
    else if (model["lattice"] == std::string("chain_lattice"))
        return new adj::ChainAdj(model["L"]);
    else if (model["lattice"] == std::string("cylinder_lattice"))
        return new adj::CylinderAdj(model["L"], model["W"]);
    else if (model["lattice"] == std::string("periodic_chain_lattice"))
        return new adj::PeriodicChainAdj(model["L"]);
    else if (model["lattice"] == std::string("periodic_ladder_lattice"))
        return new adj::PeriodicLadderAdj(model["L"]);
    else if (model["lattice"] == std::string("periodic_square_lattice"))
        return new adj::PeriodicSquareLatticeAdj(model["L"], model["W"]);
    else if (model["lattice"] == std::string("snake_square_lattice"))
        return new adj::SnakeSquareAdj(model["L"], model["W"]);
    else {
        throw std::runtime_error("Don't know this lattice!");
        return NULL;
    }
};  

template<class Matrix>
mpos::Hamiltonian<Matrix, U1> * hamil_factory(ModelParameters & model)
{
    if (model["model"] == std::string("heisenberg"))
        return new mpos::Heisenberg<Matrix>(model["Jxy"], model["Jz"]);
    else if (model["model"] == std::string("biquadratic"))
        return new mpos::Spin1BlBq<Matrix>(cos(M_PI * model["theta"]),
                                           sin(M_PI * model["theta"]),
                                           model["h0"]);
    else if (model["model"] == std::string("HCB"))
        return new mpos::HCB<Matrix>();
    else if (model["model"] == std::string("FreeFermions"))
        return new mpos::FreeFermions<Matrix>();
    else {
        throw std::runtime_error("Don't know this model!");
        return NULL;
    }
}

int main(int argc, char ** argv)
{
    Timer everything("everything");
    everything.begin();
    
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
    
    std::string chkpfile = parms["chkpfile"];
    std::string rfile = parms["resultfile"];

    bool dns = (parms["donotsave"] != 0);
    
    bool restore = false;
    {
        struct stat tmp;
        if (stat(chkpfile.c_str(), &tmp) == 0 && S_ISREG(tmp.st_mode))
        {
            maquis::cout << "Restoring state." << std::endl;
            restore = true;
        }
    }
    
    adj::Adjacency * adj = adj_factory(model);
    mpos::Hamiltonian<Matrix, grp> * H = hamil_factory<Matrix>(model);
    Index<U1> phys = H->get_phys();
    
    int total_charge = model["u1_total_charge"];
    MPS<Matrix, grp> mps(adj->size(), 5, phys, total_charge, *new default_mps_init<Matrix, grp>);
    
    int sweep = 0;
    if (restore) {
        storage::archive ar_in(chkpfile);
        ar_in["/state"] >> mps;
        ar_in["/status/sweep"] >> sweep;
        ++sweep;
    }
    
    {
        storage::archive ar(rfile, "w");
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
    }
    
    if (!dns) {
        storage::archive ar(chkpfile, "w");
        ar["/parameters"] << parms;
        ar["/parameters"] << model;
    }
    
    storage::setup(parms);
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    std::vector<double> energies, entropies;
    std::vector<std::size_t> truncations;
    
    bool early_exit = false;
    
    mpos::MPOMaker<Matrix, grp> mpom(*adj, *H);
    mpom.add_bond_ops();
    H->push_extra_terms(mpom, *adj);
    MPO<Matrix, grp> mpo = mpom.create_mpo();
    
    maquis::cout << maquis::real(expval(mps, mpo, 0)) << std::endl;
    maquis::cout << maquis::real(expval(mps, mpo, 1)) << std::endl;
    {   
        ss_optimize<Matrix, grp, storage::disk> optimizer(mps, parms);
        
        for (int sweep = 0; sweep < parms["nsweeps"]; ++sweep) {
            gettimeofday(&snow, NULL);
            
            std::pair<std::vector<double>, std::vector<std::size_t> > r;
            optimizer.sweep(mpo, sweep);
            
            entropies = calculate_bond_entropies(mps);
            
            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            
            {
                storage::archive ar(rfile, "w");
                
                std::ostringstream oss;
                oss << "/simulation/sweep" << sweep << "/results";
                ar[oss.str().c_str()] << storage::log;
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Iteration Entropies/mean/value";
                ar[oss.str().c_str()] << entropies;
                
                maquis::cout << "Sweep done after " << elapsed << " seconds." << std::endl;
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Runtime/mean/value";
                ar[oss.str().c_str()] << std::vector<double>(1, elapsed);
            }
                
            if (!dns)
            {
                storage::archive ar(chkpfile, "w");
                
                ar["/state"] << mps;
                ar["/status/sweep"] << sweep;
            }
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec); 
            int rs = parms["run_seconds"];
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
    }
    
    
    storage::disk::sync();
    
    if (!early_exit)
    {
        storage::archive ar(rfile, "w");
        
        measure(mps, *adj, *H, model, ar);
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
        maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
        
        ar["/simulation/results/Iteration Energy/mean/value"] << energies;
        ar["/simulation/results/Runtime/mean/value"] << std::vector<double>(1, elapsed);
        
        ar["/spectrum/results/Entropy/mean/value"] << entropies;
        
        everything.end();
    }
}
