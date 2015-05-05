/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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


namespace measurements_details {

    template <class symm, class = void>
    class checkpg
    {
    public:

        template <class matrix>
        bool operator()(term_descriptor<double> const & term, boost::shared_ptr<TagHandler<matrix, symm> > tag_handler) 
        {
            return true;
        }
    };

    template <class symm>
    class checkpg<symm, typename boost::enable_if<symm_traits::HasPG<symm> >::type>
    {
    public:
        typedef typename symm::charge charge;
        typedef typename symm::subcharge subcharge;

        template <class matrix>
        bool operator()(term_descriptor<double> const & term, boost::shared_ptr<TagHandler<matrix, symm> > tag_handler) 
        {
            typedef typename TagHandler<matrix, symm>::op_t op_t;

            maquis::cout <<  "number of terms..." << term.size() << std::endl; 
            op_t product = tag_handler->get_op(term.operator_tag(0));
            //maquis::cout << " product.basis().size() start --> " << product.basis().size() << std::endl;
            //maquis::cout << " product.basis().left/right() start --> " << product.basis().left_charge(0) << product.basis().right_charge(0) << std::endl;
            for (std::size_t p = 1; p < term.size(); ++p) {
                op_t tmp2 = tag_handler->get_op(term.operator_tag(p));
                //maquis::cout << " product.basis().size() current op --> " << tmp2.basis().size() << std::endl;
                //maquis::cout << " product.basis().left/right() current op  --> " << tmp2.basis().left_charge(0) << tmp2.basis().right_charge(0) << std::endl;
                op_t tmp;
                gemm(product, tmp2, tmp);
                swap(tmp, product);
                maquis::cout << "  checkpg prod " << product << std::endl;
                //maquis::cout << " product.basis().size() after gemm  --> " << product.basis().size() << std::endl;
                //if(product.basis().size() > 0)
                 //   maquis::cout << " product.basis().left/right() after gemm  --> " << product.basis().left_charge(0) << product.basis().right_charge(0) << std::endl;
            }

            maquis::cout << "checkpg prod " << product << std::endl;
            //maquis::cout << " product.basis().size() final --> " << product.basis().size() << std::endl;

            if(product.basis().size() > 0)
                for (std::size_t p = 0; p < product.basis().size(); ++p) {
                maquis::cout << " product.basis().left/right() check --> " << product.basis().left_charge(p) << product.basis().right_charge(p) << std::endl;
                    if (product.basis().left_charge(p) != product.basis().right_charge(p))
                        return false;
                maquis::cout << " done checking... returning true" << std::endl;
                }
                return true;
            return false;
        }
    };

}

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
        boost::shared_ptr<TagHandler<matrix, symm> > tag_handler = model.operators_table();

        // load state
        MPS<matrix, symm> mps;
        std::string wvf = parms["chkpfile"];
        load(wvf, mps);

        //int pos_[4] = {0, 1, 4, 5};
        //std::vector<int> pos(pos_, pos_ + 4);
        int pos_[6] = {0, 0, 2, 0, 0, 2};
        std::vector<int> pos(pos_, pos_ + 6);

        typedef typename operator_selector<matrix, symm>::type op_t;
        typedef typename OPTable<matrix, symm>::tag_type tag_type;

        //tag_type op1 = model.get_operator_tag("create_down", lattice.get_prop<symm::subcharge>("type",  pos[0]));
        //tag_type op2 = model.get_operator_tag("create_down", lattice.get_prop<symm::subcharge>("type",  pos[1]));
        //tag_type op3 = model.get_operator_tag("destroy_down", lattice.get_prop<symm::subcharge>("type", pos[2]));
        //tag_type op4 = model.get_operator_tag("destroy_down", lattice.get_prop<symm::subcharge>("type", pos[3]));
        
        //tag_type ops_[4] = {op1, op2, op3, op4};
        //std::vector<tag_type> ops(ops_, ops_ + 4);
          

          tag_type op1 = model.get_operator_tag("create_up", lattice.get_prop<symm::subcharge>("type",  pos[0]));
          tag_type op2 = model.get_operator_tag("create_down", lattice.get_prop<symm::subcharge>("type",  pos[1]));
          tag_type op3 = model.get_operator_tag("create_up", lattice.get_prop<symm::subcharge>("type",  pos[2]));
          tag_type op4 = model.get_operator_tag("destroy_up", lattice.get_prop<symm::subcharge>("type", pos[3]));
          tag_type op5 = model.get_operator_tag("destroy_down", lattice.get_prop<symm::subcharge>("type", pos[4]));
          tag_type op6 = model.get_operator_tag("destroy_up", lattice.get_prop<symm::subcharge>("type", pos[5]));

          tag_type ops_[6] = {op1, op2, op3, op4, op5, op6};
          std::vector<tag_type> ops(ops_, ops_ + 6);

        term_descriptor<double> term = generate_mpo::arrange_operators(pos, ops, tag_handler);
        maquis::cout << term << std::endl;
        for (int i=0; i<term.size(); ++i)
            maquis::cout << tag_handler->get_op(term.operator_tag(i));

        std::vector<tag_type> identities, fillings;
        identities.push_back(model.identity_matrix_tag(0));
        identities.push_back(model.identity_matrix_tag(1));
        fillings.push_back(model.filling_matrix_tag(0));
        fillings.push_back(model.filling_matrix_tag(1));

        // check if term is allowed by symmetry
        if(not measurements_details::checkpg<symm>()(term, tag_handler))
               maquis::cout << "term 0 by symmetry" << std::endl;


        MPO<matrix, symm> mpo = generate_mpo::make_1D_mpo(pos, ops, identities, fillings, tag_handler, lattice);
        double value = expval(mps, mpo);
        maquis::cout << "Expval is: " << value << std::endl; 
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
