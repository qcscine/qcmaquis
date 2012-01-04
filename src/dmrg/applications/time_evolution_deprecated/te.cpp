
#include "dmrg_version.h"

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#ifdef USE_MTM
#include "types/mt_matrix/mt_matrix.h"
#include "types/mt_matrix/algorithms.hpp"

#ifdef IMG_ONLY
typedef maquis::types::mt_matrix<double> Matrix;
#else
typedef maquis::types::mt_matrix<std::complex<double> > Matrix;
#endif

#else

#ifdef IMG_ONLY
typedef maquis::types::dense_matrix<double> Matrix;
#else
typedef maquis::types::dense_matrix<std::complex<double> > Matrix;
#endif

#endif

#include <alps/hdf5.hpp>

#include "dmrg/utils/DmrgParameters2.h"

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"
#include "dmrg/utils/random.hpp"

#include "dmrg/sim/te_utils.hpp"
#include "dmrg/mp_tensors/te.h"

#include "dmrg/models/factory.h"

using namespace app;

#ifdef UseTwoU1
typedef TwoU1 grp;
#else
typedef U1 grp;
#endif

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

std::ostream& operator<< (std::ostream& os, std::pair<grp::charge, std::size_t> const& p)
{
    os << "(" << p.first << " : " << p.second << ")";
    return os;
}

std::ostream& operator<< (std::ostream& os, index_product_iterator<grp>::value_type const& v)
{
    //std::copy(v.begin(), v.end(), std::ostream_iterator<std::pair<symm::charge, std::size_t> >(os, " "));
    for (int i=0; i<v.size(); ++i)
        os << v[i] << " ";
    return os;
}

std::ostream& operator<< (std::ostream& os, std::pair<MultiIndex<grp>::coord_t, MultiIndex<grp>::coord_t> const& p)
{
    os << p.first << ", " << p.second;
    return os;
}


template<class Matrix>
mps_initializer<Matrix, grp> * initializer_factory(BaseParameters & params)
{
    if (params.get<std::string>("init_state") == "default")
        return new default_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "const")
        return new const_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "thin")
        return new thin_mps_init<Matrix, grp>();
    else if (params.get<std::string>("init_state") == "thin_const")
        return new thin_const_mps_init<Matrix, grp>();
    else {
        throw std::runtime_error("Don't know this initial state.");
        return NULL;
    }
}

std::vector<MPO<Matrix, grp> >
getU(std::vector<Hamiltonian<Matrix, grp> > const & split_H, Lattice_ptr lat,
     double dt, bool img)
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
    
    std::vector<MPO<Matrix, grp> > expMPO(split_H.size(), MPO<Matrix, grp>(lat->size()));
    for (int i=0; i<split_H.size(); ++i)
//        expMPO[i] = make_exp_mpo(lat->size(), split_H[i], alpha);
        expMPO[i] = make_exp_mpo(lat->size(), split_H[i], alpha);
    return expMPO;
}

int main(int argc, char ** argv)
{
    cout << DMRG_VERSION_STRING << endl;
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
    dmrg_random::engine.seed(parms.get<int>("seed"));
    
#ifdef IMG_ONLY
    if (parms.get<int>("nsweeps_img") != parms.get<int>("nsweeps")) {
        cerr << "IMG_ONLY code, make sure that nsweeps_img == nsweeps." << endl;
        exit(1);
    }
#endif
    
    Lattice_ptr lat;
    model_traits<Matrix, grp>::model_ptr phys_model;
    model_parser<Matrix, grp>(parms.get<std::string>("lattice_library"), parms.get<std::string>("model_library"), model, lat, phys_model);
    Hamiltonian<Matrix, grp> H = phys_model->H();
    grp::charge initc = phys_model->initc(model);
    Measurements<Matrix, grp> measurements = phys_model->measurements();
    Index<grp> phys = H.get_phys();
    std::cout << "initc: " << initc << std::endl;
    
    std::cout << measurements << std::endl;
    
    std::vector<Hamiltonian<Matrix, grp> > split_H = separate_overlaps(H);
    std::vector<MPO<Matrix, grp> > expMPO = getU(split_H, lat, parms.get<double>("dt"), true);
    
    Measurements<Matrix, grp> meas_always;
    if (!parms.get<std::string>("always_measure").empty()) {
        meas_always.set_identity(measurements.get_identity());
        std::vector<std::string> meas_list = parms.get<std::vector<std::string> >("always_measure");
        for (int i=0; i<meas_list.size(); ++i)
            meas_always.add_term(measurements.get(meas_list[i]));
    }
    
    MPO<Matrix, grp> mpo = make_mpo(lat->size(), H);
    MPO<Matrix, grp> mpoc = mpo;
    if (parms.get<int>("use_compressed") > 0)
        mpoc.compress(1e-12);
    
    MPS<Matrix, grp> initial_mps(lat->size(),
                                 parms.get<std::size_t>("init_bond_dimension"),
                                 phys, initc,
                                 *initializer_factory<Matrix>(parms));
    
    int sweep = 0;
    if (restore) {
        alps::hdf5::archive h5ar_in(chkpfile);
        h5ar_in >> alps::make_pvp("/state", initial_mps);
        h5ar_in >> alps::make_pvp("/status/sweep", sweep);
        ++sweep;
    } else if (parms.get<std::string>("initfile").size() > 0) {
        alps::hdf5::archive h5ar_in(parms.get<std::string>("initfile"));
        h5ar_in >> alps::make_pvp("/state", initial_mps);
    }
    
    {
        alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    
    if (!dns) {
        alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
        h5ar << alps::make_pvp("/parameters", parms);
        h5ar << alps::make_pvp("/parameters", model);
        h5ar << alps::make_pvp("/version", DMRG_VERSION_STRING);
    }
    
//    StreamStorageMaster ssm(parms.get<std::string>("storagedir"));
    NoopStorageMaster ssm;
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    std::vector<double> energies, entropies, renyi2;
    std::vector<std::size_t> truncations;
    MPS<Matrix, grp> cur_mps = initial_mps, old_mps = initial_mps;
    
#ifndef MEASURE_ONLY
    
    bool early_exit = false;
    {
        for ( ; sweep < parms.get<int>("nsweeps"); ++sweep) {
            if (sweep == parms.get<int>("nsweeps_img"))
                expMPO = getU(split_H, lat, parms.get<double>("dt"), false);
            
            gettimeofday(&snow, NULL);
            
            Logger iteration_log;
            
            /* this will need some thinking about how to do it most efficiently */
            for (int which = 0; which < expMPO.size(); ++which)
            {
                old_mps = cur_mps;
                time_evolve<Matrix, grp, NoopStorageMaster> evolution(cur_mps,
                                                                      expMPO[which],
                                                                      parms, ssm);
                for (int k = 0; k < 5; ++k)
                    evolution.sweep(sweep, iteration_log);
                evolution.finalize();
                cur_mps = evolution.get_current_mps();
//                cout << "Overlap " << overlap(cur_mps, old_mps) << endl;
//                cout << norm(cur_mps) << endl;
            }
            
//            entropies = calculate_bond_entropies(cur_mps);
            
            double energy = expval(cur_mps, mpoc);
            cout << "Energy " << energy << endl;

            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);
            cout << "Sweep " << sweep << " done after " << elapsed << " seconds." << endl;
            
            {
//                alps::hdf5::archive h5ar(rfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
//                
//                std::ostringstream oss;
//                
//                oss.str("");
//                oss << "/simulation/sweep" << sweep << "/results";
//                h5ar << alps::make_pvp(oss.str().c_str(), iteration_log);
//
//                oss.str("");
//                oss << "/simulation/sweep" << sweep << "/results/Energy/mean/value";
//                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, energy));
//
//                oss.str("");
//                oss << "/simulation/sweep" << sweep << "/results/Iteration Entropies/mean/value";
//                h5ar << alps::make_pvp(oss.str().c_str(), entropies);
//                
//                oss.str("");
//                oss << "/simulation/sweep" << sweep << "/results/Runtime/mean/value";
//                h5ar << alps::make_pvp(oss.str().c_str(), std::vector<double>(1, elapsed));                
            }
            
            if (sweep % parms.get<int>("measure_each") == 0)
            {
                std::ostringstream oss;
                oss << "/simulation/sweep" << sweep << "/results/";
                if (meas_always.n_terms() > 0)
                    measure_on_mps(cur_mps, *lat, meas_always, rfile, oss.str());
            }
            
            
            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }

            if (!dns && (early_exit ||
                         sweep % parms.get<int>("chkp_each") == 0 ||
                         sweep == parms.get<int>("nsweeps")))
            {
                alps::hdf5::archive h5ar(chkpfile, alps::hdf5::archive::WRITE | alps::hdf5::archive::REPLACE);
                
                h5ar << alps::make_pvp("/state", cur_mps);
                h5ar << alps::make_pvp("/status/sweep", sweep);
            }
            
            if (early_exit)
                break;
        }
    }
#endif
        
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    cout << "Task took " << elapsed << " seconds." << endl;
}
