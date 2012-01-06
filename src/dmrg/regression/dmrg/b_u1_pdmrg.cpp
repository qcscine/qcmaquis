#ifdef MPI_PARALLEL
#include "ambient/ambient.h"
typedef ambient::dim2 dim;
#endif

#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/utils/utils.hpp"
#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/concept/matrix_interface.hpp"
#include "types/p_dense_matrix/concept/resizable_matrix_interface.hpp"
//#include "types/p_dense_matrix/p_dense_matrix_algorithms.hpp"
#include "types/p_dense_matrix/algorithms.hpp"
#include "types/p_dense_matrix/p_diagonal_matrix.h"

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#include "types/utils/matrix_cast.h"

typedef maquis::types::p_dense_matrix<double> Matrix;


#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "types/p_dense_matrix/block_matrix_callbacks.hpp"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"

#include "dmrg/mp_tensors/ss_optimize.hpp"

#include "dmrg/deprecated/mpos/alps_adjacency.h"

#include "b_measurements.h"
#include "b_adjacency.h"
#include "b_generate_mpo.h"
#include "b_hamiltonians.h"
#include "b_2u1_hamiltonians.h"
#include "b_superfermion_mpo.h"
#include "b_DmrgParameters.h"

#ifdef UseTwoU1
typedef TwoU1 grp;
#else
typedef U1 grp;
#endif

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

b_adj::Adjacency * adj_factory(BaseParameters & model)
{
    if (model.get<std::string>("lattice") == std::string("square_lattice"))
        return new b_adj::SquareAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("chain_lattice"))
        return new b_adj::ChainAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("cylinder_lattice"))
        return new b_adj::CylinderAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("periodic_chain_lattice"))
        return new b_adj::PeriodicChainAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("periodic_ladder_lattice"))
        return new b_adj::PeriodicLadderAdj(model.get<int>("L"));
    else if (model.get<std::string>("lattice") == std::string("periodic_square_lattice"))
        return new b_adj::PeriodicSquareLatticeAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("snake_square_lattice"))
        return new b_adj::SnakeSquareAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("alps_lattice"))
        return new adj::ALPSAdj(model.get<std::string>("alps_lattice"));
    else {
        throw std::runtime_error("Don't know this lattice!");
        return NULL;
    }
};  

template<class Matrix>
b_mpos::Hamiltonian<Matrix, grp> * hamil_factory(BaseParameters & model, int sweep)
{
#ifdef UseTwoU1
    if (model.get<std::string>("model") == std::string("biquadratic"))
        return new b_mpos::TwoU1_Spin1BlBq<Matrix>(model);
    else if (model.get<std::string>("model") == std::string("fermi_hubbard"))
        return new b_mpos::TwoU1_FermiHubbard<Matrix>(model.get<double>("t"), model.get<double>("U"));
    else {
        throw std::runtime_error("Don't know this model!");
        return NULL;
    }
#else
    if (model.get<std::string>("model") == std::string("heisenberg"))
        return new b_mpos::Heisenberg<Matrix>(model.get<double>("Jxy"), model.get<double>("Jz"));
    else if (model.get<std::string>("model") == std::string("biquadratic"))
        return new b_mpos::Spin1BlBq<Matrix>(model);
    else if (model.get<std::string>("model") == std::string("HCB"))
        return new b_mpos::HCB<Matrix>();
    else if (model.get<std::string>("model") == std::string("FreeFermions"))
        return new b_mpos::FreeFermions<Matrix>();
    else if (model.get<std::string>("model") == std::string("superfermion"))
        return new b_mpos::Superfermions<Matrix>(model, sweep);
    else if (model.get<std::string>("model") == std::string("alternate_superfermion"))
        return new b_mpos::AlternateSuperfermions<Matrix>(model, sweep);
    else if (model.get<std::string>("model") == std::string("AlternateFreeFermions"))
        return new b_mpos::AlternateFreeFermions<Matrix>(model);
    else {
        throw std::runtime_error("Don't know this model!");
        return NULL;
    }
#endif
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

int main(int argc, char ** argv)
{
    ambient::init();
    ambient::layout >> dim(1,1), dim(1,1), dim(1,1); 

    if (argc != 3)
    {
        zout << "Usage: <parms> <model_parms>" << endl;
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
  
    srand48(parms.get<int>("seed"));//+ambient::rank() :D
    
    b_adj::Adjacency * adj = adj_factory(model);
    b_mpos::Hamiltonian<Matrix, grp> * H = hamil_factory<Matrix>(model, 0);
    Index<grp> phys = H->get_phys();
    
    b_mpos::MPOMaker<Matrix, grp> mpom(*adj, *H);
    mpom.add_bond_ops();
    H->push_extra_terms(mpom, *adj);
    
    MPO<Matrix, grp> mpo = mpom.create_mpo();
    MPO<Matrix, grp> mpoc = mpo;
    if(parms.get<int>("use_compressed") > 0)
        mpoc.compress(1e-12);
    
#ifdef UseTwoU1
    int tc1 = model.get<int>("u1_total_charge1");
    int tc2 = model.get<int>("u1_total_charge2");
    grp::charge initc;
    initc[0] = tc1;
    initc[1] = tc2;
#else
    int initc = model.get<int>("u1_total_charge");
#endif
    MPS<Matrix, grp> mps(adj->size(),
                         parms.get<std::size_t>("init_bond_dimension"),
                         phys, initc,
                         *initializer_factory<Matrix>(parms));
    
    int sweep = 0;

/*
    std::string p = parms.get<std::string>("storagedir") + boost::lexical_cast<std::string>(ambient::rank());
*/
    //StreamStorageMaster ssm(ambient::operator+(parms.get<std::string>("storagedir"),ambient::rank()));
    NoopStorageMaster ssm;
    
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    
    std::vector<double> energies, entropies, renyi2;
    std::vector<std::size_t> truncations;

#ifndef MEASURE_ONLY
    
    bool early_exit = false;
    {   
        ss_optimize<Matrix, grp, NoopStorageMaster> optimizer(mps,
                                                              parms.get<int>("use_compressed") == 0 ? mpo : mpoc,
                                                              parms, ssm);
        
        for(; sweep < parms.get<int>("nsweeps"); ++sweep)
        {
            gettimeofday(&snow, NULL);
            Logger iteration_log;
            //assert(false);
            optimizer.sweep(sweep, iteration_log);

            zout << "\nSweep is done\n";
            
            //ssm.sync();
            
            entropies = calculate_bond_entropies(mps);
            
            gettimeofday(&sthen, NULL);
            double elapsed = sthen.tv_sec-snow.tv_sec + 1e-6 * (sthen.tv_usec-snow.tv_usec);


            gettimeofday(&then, NULL);
            elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);            
            int rs = parms.get<int>("run_seconds");
            if (rs > 0 && elapsed > rs) {
                early_exit = true;
                break;
            }
        }
        //ssm.sync();
    }
#endif
    
#ifdef MEASURE_ONLY
    if(ambient::is_master())
    {

        Timer tvn("vN entropy"), tr2("Renyi n=2");
        zout << "Calculating vN entropy." << endl;
        tvn.begin(); entropies = calculate_bond_entropies(mps); tvn.end();
        zout << "Calculating n=2 Renyi entropy." << endl;
        tr2.begin(); renyi2 = calculate_bond_renyi_entropies(mps, 2); tr2.end();
        
        double energy = expval(mps, mpoc);
        zout << "Energy before: " << expval(mps, mpo) << endl;
        zout << "Energy after: " << expval(mps, mpoc) << endl;
        h5ar << alps::make_pvp("/spectrum/results/Energy/mean/value", std::vector<double>(1, energy));
        
        if (parms.get<int>("calc_h2") > 0) {
            Timer tt1("square"), tt2("compress");
            tt1.begin(); MPO<Matrix, grp> mpo2 = square_mpo(mpo); tt1.end();
            tt2.begin(); mpo2.compress(1e-12); tt2.end();
            
            Timer t3("expval mpo2"), t4("expval mpo2c");
            
            t4.begin();
            double energy2 = expval(mps, mpo2, true);
            t4.end();
        
            zout << "Energy^2: " << energy2 << endl;
            zout << "Variance: " << energy2 - energy*energy << endl;
        }

        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
        
        zout << "Task took " << elapsed << " seconds." << endl;
    }
#endif
    zout << "Exiting... Finalize!\n";
    ambient::finalize();
    return 0;
}
