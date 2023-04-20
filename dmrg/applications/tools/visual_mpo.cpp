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


int main(int argc, char ** argv)
{
    try {
        DmrgOptions opt(argc, argv);
        if (!opt.valid) return 0;
        DmrgParameters parms = opt.parms;
        
        maquis::cout.precision(10);
        
        /// Parsing model
        Lattice lattice = Lattice(parms);
        Model<matrix, symm> model = Model<matrix, symm>(lattice, parms);
        
        MPO<matrix, symm> mpo = make_mpo(lattice, model);
        
        typedef MPOTensor<matrix, symm>::col_proxy col_proxy;
        typedef MPOTensor<matrix, symm>::index_type index_type;
        typedef OPTable<matrix, symm>::tag_type tag_type;
        
        index_type length = lattice.size();
        std::stringstream labels, edges;
        
        // as_left
        for (index_type b2=0; b2<mpo[0].col_dim(); ++b2) {
            col_proxy col_b2 = mpo[0].column(b2);
            for (col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
                index_type b1 = col_it.index();
                MPOTensor_detail::term_descriptor<matrix, symm> access = mpo[0].at(b1,b2);
                
                tag_type tag = mpo[0].tag_number(b1,b2);
                double v = access.scale;
                
                std::stringstream curr;
                curr << "M" << 0 << "_" << b1 << "_" << b2;
                
                labels << curr.str() << "[label=\"" << v << "*[" << tag << "]" << "\"]" << ";\n";
                edges << "start" << " -> " << curr.str() << "[label=" << b1 << "]" << ";\n";
            }
        }
        
        // as_bulk
        for (size_t p = 1; p < length; ++p) {
            
            for (index_type b2=0; b2<mpo[p].col_dim(); ++b2) {
                col_proxy col_b2 = mpo[p].column(b2);
                for (col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
                    index_type b1 = col_it.index();
                    MPOTensor_detail::term_descriptor<matrix, symm> access = mpo[p].at(b1,b2);
                    
                    tag_type tag = mpo[p].tag_number(b1,b2);
                    double v = access.scale;
                    
                    
                    
                    std::stringstream curr;
                    curr << "M" << p << "_" << b1 << "_" << b2;
                    
                    labels << curr.str() << "[label=\"" << v << "*[" << tag << "]" << "\"]" << ";\n";
                    
                    col_proxy col_prev = mpo[p-1].column(b1);
                    for (col_proxy::const_iterator prev_it = col_prev.begin(); prev_it != col_prev.end(); ++prev_it) {
                        index_type prev_b1 = prev_it.index();
                        
                        std::stringstream prev;
                        prev << "M" << p-1 << "_" << prev_b1 << "_" << b1;
                        
                        edges << prev.str() << " -> " << curr.str() << "[label=" << b1 << "]" << ";\n";
//                        edges << prev.str() << " -> " << curr.str() << ";\n";
                    }
                }
            }
        }
        
        
        // as_right
        for (index_type b2=0; b2<mpo[length-1].col_dim(); ++b2) {
            col_proxy col_b2 = mpo[length-1].column(b2);
            for (col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
                index_type b1 = col_it.index();
                MPOTensor_detail::term_descriptor<matrix, symm> access = mpo[length-1].at(b1,b2);
                
                tag_type tag = mpo[length-1].tag_number(b1,b2);
                double v = access.scale;
                
                std::stringstream prev;
                prev << "M" << length-1 << "_" << b1 << "_" << b2;
                
                edges << prev.str() << " -> " << "fin" << "[label=" << b2 << "]" << ";\n";
            }
        }
        
        std::cout << "digraph hamiltonian {\n";
        std::cout << labels.str() << "\n";
        std::cout << edges.str() << "\n";
        std::cout << "}\n";
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}

