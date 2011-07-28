#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/aligned_allocator.h"

#ifdef USE_GPU
#include <cublas.h>
#endif

#ifdef USE_MTM
#define USE_MTM_MAIN
#include "dense_matrix/mt_matrix.h"
typedef mt_matrix<std::complex<double> > Matrix;
#else

#ifdef IMG_ONLY
typedef blas::dense_matrix<double> Matrix;
#else
typedef blas::dense_matrix<std::complex<double> > Matrix;
#endif

#endif

#include <alps/hdf5.hpp>

#include "utils/DmrgParameters2.h"

#include "block_matrix/indexing.h"
#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"
#include "mp_tensors/contractions.h"
#include "mp_tensors/mps_mpo_ops.h"
#include "mp_tensors/mpo_ops.h"
#include "mp_tensors/mps_initializers.h"

#include "utils/stream_storage.h"
#include "utils/logger.h"

#include "te_utils.hpp"
#include "mp_tensors/evolve.h"

#include "models.h"

using namespace app;

#ifdef UseTwoU1
typedef TwoU1 grp;
#else
typedef U1 grp;
#endif

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

template<class Matrix>
mps_initializer<Matrix, grp> * initializer_factory(BaseParameters & params)
{
    if (params.get<std::string>("init_state") == "default")
        return new default_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "mott")
        return new mott_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "thin")
        return new thin_mps_init<Matrix, grp>();
    else {
        throw std::runtime_error("Don't know this initial state.");
        return NULL;
    }
}

std::vector<std::map<std::size_t, block_matrix<Matrix, grp> > >
get_U (Hamiltonian<Matrix, grp> const & H, double dt, bool img)
{
#ifdef IMG_ONLY
    double alpha = -dt;
#else
    std::complex<double> I;
    if (img)
        I = std::complex<double>(1, 0);
    else
        I = std::complex<double>(0, 1);
    std::complex<double> alpha = -I*dt;
#endif

    std::vector<Hamiltonian<Matrix, grp> > split_H = separate_overlaps(H);
    std::vector<std::map<std::size_t, block_matrix<Matrix, grp> > > expH(split_H.size());
    for (int i=0; i<split_H.size(); ++i)
        expH[i] = make_exp_nn(split_H[i], alpha);
    
    std::cout << expH.size() << " non overlapping Hamiltonians" << std::endl;
    for (int i=0; i<expH.size(); ++i)
    {
        std::cout << "Hamiltonian " << i << std::endl;
        for (std::map<std::size_t, block_matrix<Matrix, grp> >::const_iterator it = expH[i].begin();
             it != expH[i].end();
             ++it)
        {
            std::cout << " ** position " << it->first << std::endl;
            std::cout << " ** matrix " << std::endl << it->second << std::endl;
            
        }
    }
    
    return expH;
}

int main(int argc, char ** argv)
{
    if (argc != 3)
    {
        cout << "Usage: <parms> <model_parms>" << endl;
        exit(1);
    }
    
    cout.precision(10);
    
#ifdef USE_GPU
	cublasInit();
#endif
    
    std::ifstream param_file(argv[1]);
    if (!param_file) {
        cerr << "Could not open parameter file." << endl;
        exit(1);
    }
    DmrgParameters parms(param_file);
    
    std::string model_file(argv[2]);
    
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
    
    srand48(parms.get<int>("seed"));
    
#ifdef IMG_ONLY
    if (parms.get<int>("nsweeps_img") != parms.get<int>("nsweeps")) {
        cerr << "IMG_ONLY code, make sure that nsweeps_img == nsweeps." << endl;
        exit(1);
    }
#endif
    
    Lattice * lat;
    Hamiltonian<Matrix, grp> H;
    grp::charge initc;
    Measurements<Matrix, grp> measurements;
    ModelParameters model = model_parser(parms.get<std::string>("model_library"), model_file, lat, H, initc, measurements);
    Index<grp> phys = H.get_phys();
    std::cout << "initc: " << initc << std::endl;
    
//    std::cout << "Hamiltonian: " << std::endl << H << std::endl;
//    std::cout << "Measurements: " << std::endl << measurements << std::endl;
        
    Measurements<Matrix, grp> meas_always;
    if (!parms.get<std::string>("always_measure").empty()) {
        meas_always.set_identity(measurements.get_identity());
        std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
        for (int i=0; i<meas_list.size(); ++i)
            meas_always.add_term(measurements.get(meas_list[i]));
    }

    MPO<Matrix, grp> mpo = make_mpo(lat->size(), H);
    
    MPS<Matrix, grp> mps(lat->size(),
                         parms.get<std::size_t>("init_bond_dimension"),
                         phys, initc,
                         *initializer_factory<Matrix>(parms));
    
    int sweep = 0;
    if (restore) {
        alps::hdf5::archive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", mps);
        h5ar_in >> alps::make_pvp("/status/sweep", sweep);
        ++sweep;
    } else if (parms.get<std::string>("initfile").size() > 0) {
        alps::hdf5::archive h5ar_in(parms.get<std::string>("initfile"));
        h5ar_in >> alps::make_pvp("/state", mps);
    }
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
    if (!dns) {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
    }
    
//    StreamStorageMaster ssm(parms.get<std::string>("storagedir"));
    NoopStorageMaster ssm;
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    std::vector<double> energies, entropies, renyi2;
    std::vector<std::size_t> truncations;
    
    
    std::vector<std::map<std::size_t, block_matrix<Matrix, grp> > > expH;
    if (sweep < parms.get<int>("nsweeps_img"))
        expH = get_U(H, parms.get<double>("dt"), true);
    else
        expH = get_U(H, parms.get<double>("dt"), false);
    
    bool early_exit = false;
    {   
        
        for ( ; sweep < parms.get<int>("nsweeps"); ++sweep) {
            gettimeofday(&snow, NULL);
            
            Logger iteration_log;
            
            // switching from ground state to time evolution
            if (sweep == parms.get<int>("nsweeps_img"))
                expH = get_U(H, parms.get<double>("dt"), false);
            
            for (int i=0; i < expH.size(); ++i)
                mps = evolve(mps, expH[i],
                             parms.get<std::size_t>("max_bond_dimension"), parms.get<double>("truncation_final"));
            
//            entropies = calculate_bond_entropies(mps);
            
//            std::complex<double> mps_norm = norm(mps);
//            std::cout << "Norm " << mps_norm << std::endl;

            double energy = expval(mps, mpo);
            std::cout << "Energy " << energy << std::endl;
            
            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            
            if (sweep % parms.get<int>("measure_each") == 0)
            {
                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE);
                
                std::ostringstream oss;
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results";
                h5ar << alps::make_pvp(oss.str().c_str(), iteration_log);
                
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Energy/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, energy));

                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Iteration Entropies/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), entropies);
                
                cout << "Sweep " << sweep << " done after " << elapsed << " seconds." << endl;
                oss.str("");
                oss << "/simulation/sweep" << sweep << "/results/Runtime/mean/value";
                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, elapsed));                
            }
            
            if (sweep % parms.get<int>("measure_each") == 0)
            {
                std::ostringstream oss;
                oss << "/simulation/sweep" << sweep << "/results/";
                if (meas_always.n_terms() > 0)
                    measure(mps, *lat, meas_always, rfile, oss.str());
            }
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
        ssm.sync();
    }
    
    if (!dns)
    {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE);
        
        h5ar << alps::make_pvp("/state", mps);
        h5ar << alps::make_pvp("/status/sweep", sweep);
    }
    
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    cout << "Task took " << elapsed << " seconds." << endl;
    
#ifdef USE_GPU
	cublasShutdown();
#endif
}
