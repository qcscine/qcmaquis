/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/
#ifdef USE_AMBIENT
#include <mpi.h>
#endif
#include <cmath>
#include <iterator>
#include <iostream>
#include <fstream>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#ifdef USE_AMBIENT
#include "dmrg/block_matrix/detail/ambient.hpp"
typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double> matrix;
#endif

#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#if defined(USE_TWOU1)
typedef TwoU1 symm;
#elif defined(USE_NONE)
typedef TrivialGroup symm;
#elif defined(USE_U1)
typedef U1 symm;
#endif

#include "dmrg/utils/DmrgParameters2.h"


int main(int argc, char ** argv)
{
    try {
        if (argc != 3)
        {
            maquis::cout << "Usage: <parms> <model_parms>" << std::endl;
            exit(1);
        }
        
        maquis::cout.precision(10);
        
        /// Loading parameters
        std::ifstream param_file(argv[1]);
        if (!param_file) {
            maquis::cerr << "Could not open parameter file." << std::endl;
            exit(1);
        }
        DmrgParameters parms(param_file);
        
        /// Loading model
        std::ifstream model_file(argv[2]);
        if (!model_file) {
            maquis::cerr << "Could not open model file." << std::endl;
            exit(1);
        }
        ModelParameters model_parms(model_file);
        
        
        /// Parsing model
        Lattice lattice = Lattice(parms, model_parms);
        Model<matrix, symm> model = Model<matrix, symm>(lattice, parms, model_parms);
        
        MPO<matrix, symm> mpo = make_mpo(lattice, model, model_parms);
        
        for (int p = 0; p < lattice.size(); ++p) {
            std::ofstream ofs(std::string("mpo_stats."+boost::lexical_cast<std::string>(p)+".dat").c_str());
            for (int b1 = 0; b1 < mpo[p].row_dim(); ++b1) {
                for (int b2 = 0; b2 < mpo[p].col_dim(); ++b2) {
                    if (mpo[p].has(b1, b2)) ofs << mpo[p].tag_number(b1,b2)+1 << " ";
                    else ofs << "0 ";
                }
                ofs << std::endl;
            }
        }
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
