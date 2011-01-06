#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
typedef blas::dense_matrix<double> Matrix;

#include <alps/hdf5.hpp>

#include "block_matrix/indexing.h"
#include "mp_tensors/mps.h"
#include "mp_tensors/mpo.h"
#include "mp_tensors/contractions.h"
#include "mp_tensors/mps_mpo_ops.h"

#include "mp_tensors/special_mpos.h"

#include "mp_tensors/ss_optimize.h"

#include "mpos/adjancency.h"
#include "mpos/generate_mpo.h"

#include "utils/DmrgParameters.h"

typedef U1 grp;

typedef std::vector<MPOTensor<Matrix, grp> > mpo_t;
typedef Boundary<Matrix, grp> boundary_t;

Adjacency * adj_factory(ModelParameters & model)
{
    if (model.get<std::string>("lattice") == std::string("square_lattice"))
        return new SquareAdj(model.get<int>("L"), model.get<int>("W"));
    else if (model.get<std::string>("lattice") == std::string("chain_lattice"))
        return new ChainAdj(model.get<int>("L"));
    else
        return NULL;
};  

template<class Matrix>
Hamiltonian<Matrix, U1> * hamil_factory(ModelParameters & model)
{
    if (model.get<std::string>("model") == std::string("heisenberg"))
        return new Heisenberg<Matrix>(model.get<double>("Jxy"), model.get<double>("Jz"));
    else if (model.get<std::string>("model") == std::string("biquadratic"))
        return new Spin1BlBq<Matrix>(cos(M_PI * model.get<double>("theta")),
                                     sin(M_PI * model.get<double>("theta")));
    else
        return NULL;
}

int main(int argc, char ** argv)
{
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
    
    alps::hdf5::oarchive h5ar(argv[3]);
    h5ar << alps::make_pvp("/parameters", parms);
    h5ar << alps::make_pvp("/parameters", model);
    
    Adjacency * adj = adj_factory(model);
    Hamiltonian<Matrix, grp> * H = hamil_factory<Matrix>(model);
    Index<U1> phys = H->get_phys();
    
    MPS<Matrix, grp> mps(adj->size(), 5, phys);
    MPO<Matrix, grp> mpo = mpos::MPOMaker<Matrix, grp>::create_mpo(*adj, *H);
    
    cout << expval(mps, mpo, 0) << endl;
    cout << expval(mps, mpo, 1) << endl;
    
    timeval now, then;
    gettimeofday(&now, NULL);
    std::vector<double> energies = ss_optimize<Matrix, grp>(mps, mpo, parms);
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    
    cout << "Task took " << elapsed << " seconds." << endl;
    
    h5ar << alps::make_pvp("/simulation/results/Iteration Energy/mean/value", energies);
    h5ar << alps::make_pvp("/simulation/results/Runtime/mean/value", std::vector<double>(1, elapsed));
    h5ar << alps::make_pvp("/spectrum/results/Energy/mean/value", std::vector<double>(1, *energies.rbegin()));
    
}
