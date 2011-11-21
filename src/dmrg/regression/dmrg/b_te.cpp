#include <complex>
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
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"
//typedef maquis::types::dense_matrix<double, std::vector<double, aligned_allocator<double> > > Matrix;
typedef maquis::types::dense_matrix<std::complex<double> > Matrix;

#include <alps/hdf5.hpp>
//template<class T, class A>
//void save(alps::hdf5::archive & ar,
//          std::string const & p,
//          std::vector<T, A> const & v,
//          std::vector<std::size_t> size = std::vector<std::size_t>(),
//          std::vector<std::size_t> chunk = std::vector<std::size_t>(),
//          std::vector<std::size_t> offset = std::vector<std::size_t>())
//{
//    std::vector<T> foo(v.begin(), v.end());
//    ar << alps::make_pvp(p, foo);
//}
//
//template<class T, class A>
//void load(alps::hdf5::archive & ar,
//          std::string const & p,
//          std::vector<T, A> & v,
//          std::vector<std::size_t> size = std::vector<std::size_t>(),
//          std::vector<std::size_t> chunk = std::vector<std::size_t>(),
//          std::vector<std::size_t> offset = std::vector<std::size_t>())
//{
//    std::vector<T> foo;
//    ar >> alps::make_pvp(p, foo);
//    v.resize(foo.size());
//    std::copy(foo.begin(), foo.end(), v.begin());
//}

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/utils/stream_storage.h"
#include "dmrg/utils/logger.h"

#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/evolve.h"

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

block_matrix<Matrix, grp> get_hb(double dt, bool img)
{
    std::complex<double> i;
    if (img)
        i = std::complex<double>(1, 0);
    else
        i = std::complex<double>(0, 1);
    
    block_matrix<Matrix, grp> ret;
    Matrix m2(1, 1), m0(2, 2);
    
    m2(0,0) = exp(-i*dt/4.0);
    ret.insert_block(m2, -2, -2);
    ret.insert_block(m2, 2, 2);
    
    m0(1, 1) = m0(0, 0) = 0.5*exp(-i*dt/4.0)*(1.0+exp(i*dt));
    m0(0, 1) = m0(1, 0) = -0.5*exp(-i*dt/4.0)*(exp(i*dt)-1.0);
    ret.insert_block(m0, 0, 0);
    
    cout << ret << endl;
    
    return ret;
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
    b_DmrgParameters parms(param_file);
    
    std::ifstream model_file(argv[2]);
    if (!model_file) {
        cerr << "Could not open model file." << endl;
        exit(1);
    }
    b_ModelParameters model(model_file);
    
    std::string chkpfile = parms.get<std::string>("chkpfile");
    std::string rfile = parms.get<std::string>("resultfile");
    
    srand48(parms.get<int>("seed"));
    
    b_adj::Adjacency * adj = new b_adj::ChainAdj(model.get<int>("L"));
    b_mpos::Hamiltonian<Matrix, grp> * H = new b_mpos::Heisenberg<Matrix>(model.get<double>("Jxy")/4,
                                                                          model.get<double>("Jz")/4);
    Index<grp> phys = H->get_phys();
    
    b_mpos::MPOMaker<Matrix, grp> mpom(*adj, *H);
    mpom.add_bond_ops();
    H->push_extra_terms(mpom, *adj);
    MPO<Matrix, grp> mpo = mpom.create_mpo();
    
    int initc = model.get<int>("u1_total_charge");
    MPS<Matrix, grp> mps(adj->size(),
                         parms.get<std::size_t>("init_bond_dimension"),
                         phys, initc,
                         *initializer_factory<Matrix>(parms));
    
    cout << norm(mps) << endl;
    
    double dt = 0.01;
    block_matrix<Matrix, grp> op = get_hb(dt, true);
    
    for (int i = 0; i < 10000; ++i)
    {
        if (i == 100)
            op = get_hb(dt, false);
        cout << expval(mps, mpo) << endl;
        mps = evolve(mps, op, 100, 1e-10);
//        cout << norm(mps) << endl;
    }
    
//    mps = compression::l2r_compress(mps, 100, 1e-8);
//    cout << norm(mps) << endl;
}
