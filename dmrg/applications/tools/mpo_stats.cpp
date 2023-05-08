/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
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

#include "dmrg/sim/matrix_types.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/ts_ops.h"
#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"

#if defined(USE_TWOU1)
typedef TwoU1 symm;
#elif defined(USE_U1DG)
typedef U1DG symm;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG symm;
#elif defined(USE_SU2U1)
typedef SU2U1 symm;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG symm;
#elif defined(USE_NONE)
typedef TrivialGroup symm;
#elif defined(USE_U1)
typedef U1 symm;
#endif

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

std::string operator * (std::string s, int m)
{
    std::string ret("");
    for (int i=0; i < m; ++i) ret += s;
    return ret;
}

template <class Matrix, class SymmGroup>
void write_mpo(MPO<Matrix, SymmGroup> const & mpo, std::string filename, bool save_space) 
{        
    std::string space(" ");

    for (int p = 0; p < mpo.size(); ++p) {
        std::ofstream ofs(std::string(filename+boost::lexical_cast<std::string>(p)+".dat").c_str());

        typename MPOTensor<Matrix, SymmGroup>::op_table_ptr op_table = mpo[p].get_operator_table();
        unsigned maxtag = op_table->size();
        int padding = 2;
        if (maxtag < 100 || save_space) padding = 1;

        for (int b1 = 0; b1 < mpo[p].row_dim(); ++b1) {
            for (int b2 = 0; b2 < mpo[p].col_dim(); ++b2) {
                if (mpo[p].has(b1, b2))
                {
                    MPOTensor_detail::term_descriptor<Matrix, SymmGroup, true> access = mpo[p].at(b1,b2);
                    int tag = mpo[p].tag_number(b1, b2, 0);
                    if (access.size() > 1)
                        ofs << space*(padding-1) << "X" << access.size();
                    else if (tag < 10)
                        ofs << space*padding << tag;
                    else if (tag < 100)
                        ofs << space*(padding-1) << tag;
                    else
                        if (save_space)
                            if (tag%100 < 10)
                                ofs << space*padding << tag%100;
                            else
                                ofs << tag%100;
                        else
                            ofs << tag;
                }
                else ofs << space*padding << ".";
            }
            ofs << std::endl;
        }
        
        ofs << std::endl;
        
        for (unsigned tag=0; tag<op_table->size(); ++tag) {
            ofs << "TAG " << tag << std::endl;
            ofs << " * op :\n" << (*op_table)[tag] << std::endl;
        }
    }
}

int main(int argc, char ** argv)
{
    try {
        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;
        
        maquis::cout.precision(10);

        bool save_space = true;
        if (parms.defined("save_space") && !parms["save_space"])
            save_space = false;
        
        /// Parsing model
        Lattice lattice = Lattice(parms);
        Model<matrix, symm> model = Model<matrix, symm>(lattice, parms);
        
        MPO<matrix, symm> mpo = make_mpo(lattice, model);
        write_mpo(mpo, "mpo_stats.", save_space);

        MPS<matrix, symm> mps = MPS<matrix, symm>(lattice.size(), *(model.initializer(lattice, parms)));
        MPO<matrix, symm> ts_mpo;
        make_ts_cache_mpo(mpo, ts_mpo, mps);

        write_mpo(ts_mpo, "ts_mpo_stats.", save_space);
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
