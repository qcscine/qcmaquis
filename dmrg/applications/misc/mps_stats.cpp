#include <cmath>
#include <iterator>
#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/block_matrix/detail/alps.hpp"

typedef alps::numeric::matrix<double> Matrix;

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
            std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
            return 1;
        }
        MPS<Matrix, grp> mps;
        
        {
            storage::archive ar(argv[1]);
            ar["/state"] >> mps;
        }
        
        for (int i=0; i<mps.length(); ++i) {
            std::string fname = "mps_stats."+boost::lexical_cast<std::string>(i)+".dat";
            std::ofstream ofs(fname.c_str());
            mps[i].make_left_paired();
            for (int k=0; k<mps[i].data().n_blocks(); ++k)
                ofs << num_rows(mps[i].data()[k]) << "    " << num_rows(mps[i].data()[k]) << std::endl;
        }
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
