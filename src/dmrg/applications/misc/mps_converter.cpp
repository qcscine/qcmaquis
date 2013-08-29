#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/block_matrix/detail/ambient.hpp"

typedef alps::numeric::matrix<double> alps_matrix;
typedef ambient::numeric::matrix<double> ambient_matrix;

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#ifdef USE_TWOU1
typedef TwoU1 grp;
#else
#ifdef USE_NONE
typedef TrivialGroup grp;
#else
typedef U1 grp;
#endif
#endif


int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " <in.h5> <out.h5>" << std::endl;
            return 1;
        }
        MPS<ambient_atrix, grp> mps_in;
        
        {
            storage::archive ar(argv[1]);
            ar["/state"] >> mps_in;
        }
        
        MPS<alps_matrix, grp> mps_out(mps_in.length());
        for (int i=0; i<mps_in.length(); ++i) {
            mps_in[i].make_left_paired();
            block_matrix<alps_matrix, grp> m;
            for (int k=0; k<mps[i].data().n_blocks(); ++k)
                m.insert_block(maquis::bindings::convert<alps_matrix>(mps[i].data()[k]), mps[i].data().left_basis()[k], mps[i].data().right_basis()[k]);
            mps_out[i] = MPSTensor<alps_matrix, grp>(mps_in[i].site_dim(),mps_in[i].row_dim(), mps_in[i].col_dim(), m, LeftPaired);
        }
        
        {
            storage::archive ar(argv[2]);
            ar["/state"] << mps_out;
        }
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
