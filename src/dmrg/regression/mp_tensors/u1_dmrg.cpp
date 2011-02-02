#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "utils/zout.hpp"

#ifdef MPI_PARALLEL

#include "ambient/ambient.h"

#include "p_dense_matrix/p_dense_matrix.h"
#include "p_dense_matrix/concept/matrix_interface.hpp"
#include "p_dense_matrix/concept/resizable_matrix_interface.hpp"
#include "p_dense_matrix/p_dense_matrix_algorithms.h"
typedef blas::p_dense_matrix<double> Matrix;

#else

#include "dense_matrix/aligned_allocator.h"
#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"

typedef blas::dense_matrix<double > Matrix;
//typedef blas::dense_matrix<double, std::vector<double, aligned_allocator<double> > > Matrix;

#endif

#include <alps/hdf5.hpp>

#include "block_matrix/indexing.h"
#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"
#include "mp_tensors/contractions.h"
#include "mp_tensors/mps_mpo_ops.h"

#include "mp_tensors/ss_optimize.h"
#include "mpos/measurements.h"

#include "mpos/adjancency.h"
#include "mpos/generate_mpo.h"
#include "mpos/hamiltonians.h"

#include "utils/DmrgParameters.h"
#include "utils/timings.h"

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

Adjacency * adj_factory(ModelParameters & model)
{
    if (model.get<std::string>("lattice") == std::string("square_lattice"))
        return new SquareAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("chain_lattice"))
        return new ChainAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("cylinder_lattice"))
        return new CylinderAdj(model.get<int>("L"), model.get<int>("W"));
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
    
#ifdef MPI_PARALLEL
    ambient::instance().init();
#endif
    
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
    
    Adjacency * adj = adj_factory(model);
    mpos::Hamiltonian<Matrix, grp> * H = hamil_factory<Matrix>(model);
    Index<U1> phys = H->get_phys();
    
    int total_charge = model.get<int>("u1_total_charge");
    MPS<Matrix, grp> mps(adj->size(), 5, phys, total_charge);
    
    int sweep = 0;
    if (restore) {
        alps::hdf5::iarchive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", mps);
        h5ar_in >> alps::make_pvp("/status/sweep", sweep);
        ++sweep;
    }
    
    {
        alps::hdf5::oarchive h5ar(rfile);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
    if (!dns) {
        alps::hdf5::oarchive h5ar(chkpfile);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
    std::vector<double> energies, entropies;
    
    bool early_exit = false;
    
    mpos::MPOMaker<Matrix, grp> mpom(*adj, *H);
    mpom.add_bond_ops();
    H->push_extra_terms(mpom, *adj);
    MPO<Matrix, grp> mpo = mpom.create_mpo();
    
    zout << expval(mps, mpo, 0) << endl;
    zout << expval(mps, mpo, 1) << endl;
    
    timeval now, then;
    gettimeofday(&now, NULL);
    {   
        ss_optimize<Matrix, grp> optimizer(mps, parms);
        
        for (int sweep = 0; sweep < parms.get<int>("nsweeps"); ++sweep) {
            energies = optimizer.sweep(mpo, sweep);
            entropies = calculate_bond_entropies(mps);
            
            gettimeofday(&then, NULL);
            double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
            
#ifdef MPI_PARALLEL
            if(ambient::scheduler::instance().is_ambient_master()) {
#endif
                {
                    alps::hdf5::oarchive h5ar(rfile);
                    
                    std::ostringstream oss;
                    oss << "/simulation/results/sweep" << sweep << "/Iteration Energy/mean/value";
                    h5ar << alps::make_pvp(oss.str().c_str(), energies);
                    
                    oss.str("");
                    oss << "/simulation/results/sweep" << sweep << "/Iteration Entropies/mean/value";
                    h5ar << alps::make_pvp(oss.str().c_str(), entropies);
                    
                    cout << "Sweep done after " << elapsed << " seconds." << endl;
                    oss.str("");
                    oss << "/simulation/results/sweep" << sweep << "/Runtime/mean/value";
                    h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, elapsed));
                }
                
                if (parms.get<int>("donotsave") == 0)
                {
                    alps::hdf5::oarchive h5ar(chkpfile);
                    
                    h5ar << alps::make_pvp("/state", mps);
                    h5ar << alps::make_pvp("/status/sweep", sweep);
                }
                
#ifdef MPI_PARALLEL
            }
#endif
            
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
    }
    
    if (!early_exit)
    {
        alps::hdf5::oarchive h5ar(rfile);
        
        measure(mps, *adj, *H, model, h5ar);
        
        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
        zout << "Task took " << elapsed << " seconds." << endl;
        
#ifdef MPI_PARALLEL
        if(ambient::scheduler::instance().is_ambient_master()){
#endif
            h5ar << alps::make_pvp("/simulation/results/Iteration Energy/mean/value", energies);
            h5ar << alps::make_pvp("/simulation/results/Runtime/mean/value", std::vector<double>(1, elapsed));
            
            h5ar << alps::make_pvp("/spectrum/results/Entropy/mean/value", entropies);
            h5ar << alps::make_pvp("/spectrum/results/Energy/mean/value", std::vector<double>(1, *energies.rbegin()));
#ifdef MPI_PARALLEL
        }
#endif
        
#ifdef MPI_PARALLEL
        ambient::instance().finalize();
#endif
        
        everything.end();
    }
}