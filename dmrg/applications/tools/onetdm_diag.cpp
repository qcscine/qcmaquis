/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Physical Chemistry, ETH Zurich
 *               2015-2015 by Stefan Knecht <stefan.knecht@phys.chem.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
 * available from http://alps.comp-phys.org/.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
 * DEALINGS IN THE SOFTWARE.
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

#include <alps/hdf5.hpp>

#include "dmrg/block_matrix/indexing.h"
#include "dmrg/block_matrix/block_matrix.h"

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/models/model.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"

#if defined(USE_TWOU1)
typedef TwoU1 symm;
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

struct corr_measurement {
    typedef operator_selector<matrix, symm>::type op_t;
    
    std::string name;
    std::vector<op_t> op1;
    bool fermionic;
};


int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " inputfile " << std::endl;
            return 1;
        }

        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;

        std::cout << "parameter " << parms << std::endl;

        typedef operator_selector<matrix, symm>::type op_t;
        typedef Model<matrix, symm>::table_ptr table_ptr;
        typedef Model<matrix, symm>::tag_type tag_type;
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        typedef typename matrix::value_type value_type;
        value_type result;
        std::vector<typename MPS<matrix, symm>::scalar_type> vector_results;
        std::vector<std::string> labels;

        
        /// Parsing model
        Lattice lattice = Lattice(parms);
        Model<matrix, symm> model = Model<matrix, symm>(lattice, parms);
        int ntypes = lattice.maximum_vertex_type()+1;
        table_ptr tag_handler = model.operators_table();

        // types are point group irreps...
        std::cout << "ntypes " << ntypes << std::endl;

        maquis::cout.precision(15);
        std::size_t L = lattice.size();

        // load bra (chkp_lhs) and ket (chkpfile) MPS
        MPS<matrix, symm> mps1, mps2;
        load(parms["chkp_lhs"].str(),mps2);
        load(parms["chkpfile"].str(),mps1);

        /// identities and fillings
        std::vector<op_t> identities(ntypes), fillings(ntypes);
        for (size_t type = 0; type < ntypes; ++type) {
            identities[type] = model.identity_matrix(type);
            fillings[type]   = model.filling_matrix(type);
        }

        // we need to measure for local ops: c+up cup, c+down cdown, c+up cdown, c+down cup
        for (std::size_t p = 0; p < L; ++p)
        {

        vector_results.clear();
        labels.clear();
        vector_results.reserve(vector_results.size() + L);
        labels.reserve(labels.size() + L);
        {
            for (int i = 0; i < 4; ++i){
                    op_vec meas_op;
                    //if (i == 1)
                        //meas_op = count_up_op;
                    //else if (i == 1)
                        //meas_op = count_down_op;
                    //else if (i == 2)
                        //meas_op = swap_u2d_op;
                    //else if (i == 3)
                        //meas_op = swap_d2u_op;
                    //else
                     //   throw std::runtime_error("Invalid observable\n");
            }
        }

              int type = lattice.get_prop<int>("type", p);
              std::cout << "type for site p " << p << ": "<< lattice.get_prop<int>("type", p)  << std::endl;
              std::cout << "tag type for nup at site p " << p << ": "<< model.get_operator_tag("count_up", type)  << std::endl;
              std::cout << "tag type for ndown at site p " << p << ": "<< model.get_operator_tag("count_down", type)  << std::endl;
              std::cout << "tag type for u2d at site p " << p << ": "<< model.get_operator_tag("u2d", type)  << std::endl;
              std::cout << "tag type for d2u at site p " << p << ": "<< model.get_operator_tag("d2u", type)  << std::endl;

              // generate MPOs
              //std::map<std::string, MPO<matrix, symm> > mpos;

              //vector_results.push_back(val);
              std::cout << "label for site p " << p << ": "<< lattice.get_prop<std::string>("label", p)  << std::endl;

              //double val = expval(mps1, mps2, mpo)
              labels.push_back( lattice.get_prop<std::string>("label", p) );

              // save the data...
              {
                  alps::hdf5::archive oh5(parms["resultfile"].str(), "w");
                  oh5["/spectrum/results/1-TDM diagonal/mean/value"] << vector_results;
                  oh5["/spectrum/results/1-TDM diagonal/labels"] << labels;
            }
        }
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << e.what() << std::endl;
        return 1;
    }
}
