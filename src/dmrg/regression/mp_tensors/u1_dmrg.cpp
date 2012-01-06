#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "utils/zout.hpp"


#include "types/dense_matrix/aligned_allocator.h"
#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"

typedef maquis::types::dense_matrix<double > Matrix;
//typedef maquis::types::dense_matrix<double, std::vector<double, aligned_allocator<double> > > Matrix;


#include <alps/hdf5.hpp>

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

#include "dmrg/utils/stream_storage.h"

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

adj::Adjacency * adj_factory(ModelParameters & model)
{
    if (model.get<std::string>("lattice") == std::string("square_lattice"))
        return new adj::SquareAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("chain_lattice"))
        return new adj::ChainAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("cylinder_lattice"))
        return new adj::CylinderAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("periodic_chain_lattice"))
        return new adj::PeriodicChainAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("periodic_ladder_lattice"))
        return new adj::PeriodicLadderAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("periodic_square_lattice"))
        return new adj::PeriodicSquareLatticeAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("snake_square_lattice"))
        return new adj::SnakeSquareAdj(model.get<int>("L"), model.get<int>("W"));
    else {
        throw std::runtime_error("Don't know this lattice!");
        return NULL;
    }
};  

template<class Matrix>
mpos::Hamiltonian<Matrix, U1> * hamil_factory(ModelParameters & model)
{
    if (model.get<std::string>("model") == std::string("heisenberg"))
        return new mpos::Heisenberg<Matrix>(model.get<double>("Jxy"), model.get<double>("Jz"));
    else if (model.get<std::string>("model") == std::string("biquadratic"))
        return new mpos::Spin1BlBq<Matrix>(cos(M_PI * model.get<double>("theta")),
                                           sin(M_PI * model.get<double>("theta")),
                                           model.get<double>("h0"));
    else if (model.get<std::string>("model") == std::string("HCB"))
        return new mpos::HCB<Matrix>();
    else if (model.get<std::string>("model") == std::string("FreeFermions"))
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
        cout << "Usage: <parms> <model_parms>" << endl;
        exit(1);
    }
    
    zout.precision(10);
    
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
    
    std::string chkpfile = parms.get<std::string>("chkpfile");
    std::string rfile = parms.get<std::string>("resultfile");
    bool dns = (parms.get<int>("donotsave") != 0);
    
    bool restore = false;
    {
        struct stat tmp;
        if (stat(chkpfile.c_str(), &tmp) == 0 && S_ISREG(tmp.st_mode))
        {
            cout << "Restoring state." << endl;
            restore = true;
        }
    }
    
    adj::Adjacency * adj = adj_factory(model);
    mpos::Hamiltonian<Matrix, grp> * H = hamil_factory<Matrix>(model);
    Index<U1> phys = H->get_phys();
    
    int total_charge = model.get<int>("u1_total_charge");
    MPS<Matrix, grp> mps(adj->size(), 5, phys, total_charge, *new default_mps_init<Matrix, grp>);
    
    int sweep = 0;
    if (restore) {
        alps::hdf5::archive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", mps);
        h5ar_in >> alps::make_pvp("/status/sweep", sweep);
        ++sweep;
    }
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
    if (!dns) {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
    StreamStorageMaster ssm(parms.get<std::string>("storagedir").c_str());
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    std::vector<double> energies, entropies;
    std::vector<std::size_t> truncations;
    
    bool early_exit = false;
    
    mpos::MPOMaker<Matrix, grp> mpom(*adj, *H);
    mpom.add_bond_ops();
    H->push_extra_terms(mpom, *adj);
    MPO<Matrix, grp> mpo = mpom.create_mpo();
    
    zout << expval(mps, mpo, 0) << endl;
    zout << expval(mps, mpo, 1) << endl;
    
    {   
        ss_optimize<Matrix, grp, StreamStorageMaster> optimizer(mps, parms, ssm);
        
        for (int sweep = 0; sweep < parms.get<int>("nsweeps"); ++sweep) {
            gettimeofday(&snow, NULL);
            
            Logger iteration_log;
            
            std::pair<std::vector<double>, std::vector<std::size_t> > r;
            optimizer.sweep(mpo, sweep, iteration_log);
            
            entropies = calculate_bond_entropies(mps);
            
            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                std::ostringstream oss;
                oss << "/simulation/sweep" << sweep << "/results";
                h5ar << alps::make_pvp(oss.str().c_str(), iteration_log);
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Iteration Entropies/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), entropies);
                
                cout << "Sweep done after " << elapsed << " seconds." << endl;
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Runtime/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, elapsed));
            }
                
                if (parms.get<int>("donotsave") == 0)
                {
                    alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                    
                    h5ar << alps::make_pvp("/state", mps);
                    h5ar << alps::make_pvp("/status/sweep", sweep);
                }
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec); 
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
    }
    
    
    ssm.sync();
    
    if (!early_exit)
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        
        measure(mps, *adj, *H, model, h5ar);
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
        zout << "Task took " << elapsed << " seconds." << endl;
        
        h5ar << alps::make_pvp("/simulation/results/Iteration Energy/mean/value", energies);
        h5ar << alps::make_pvp("/simulation/results/Runtime/mean/value", std::vector<double>(1, elapsed));
        
        h5ar << alps::make_pvp("/spectrum/results/Entropy/mean/value", entropies);
        
        everything.end();
    }
}
