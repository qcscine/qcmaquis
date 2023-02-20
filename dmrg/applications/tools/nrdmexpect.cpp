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

MPO<matrix, symm> rdm2(Model<matrix, symm> const & model, Lattice const & lattice)
{
        typedef operator_selector<matrix, symm>::type op_t;
        typedef OPTable<matrix, symm>::tag_type tag_type;
        typedef matrix::value_type value_type;

        std::shared_ptr<TagHandler<matrix, symm> > tag_handler = model.operators_table();

        int pos_[4] = {0, 1, 2, 3};
        std::vector<int> pos(pos_, pos_ + 4);

        // B
        tag_type op1 = model.get_operator_tag("create_fill", lattice.get_prop<symm::subcharge>("type", pos[0]));
        tag_type op2 = model.get_operator_tag("create", lattice.get_prop<symm::subcharge>("type", pos[1]));
        tag_type op3 = model.get_operator_tag("destroy_fill", lattice.get_prop<symm::subcharge>("type", pos[2]));
        tag_type op4 = model.get_operator_tag("destroy", lattice.get_prop<symm::subcharge>("type", pos[3]));

        tag_type ops_[4] = {op1, op2, op3, op4};
        std::vector<tag_type> ops(ops_, ops_ + 4);

        term_descriptor<value_type> term = generate_mpo::arrange_operators(pos, ops, tag_handler);
        maquis::cout << term << std::endl;
        for (int i=0; i<term.size(); ++i)
            maquis::cout << tag_handler->get_op(term.operator_tag(i));

        std::vector<tag_type> identities, fillings;
        identities.push_back(model.identity_matrix_tag(0));
        identities.push_back(model.identity_matrix_tag(1));
        fillings.push_back(model.filling_matrix_tag(0));
        fillings.push_back(model.filling_matrix_tag(1));

        MPO<matrix, symm> mpo = generate_mpo::make_1D_mpo(pos, ops, identities, fillings, tag_handler, lattice);

        return mpo;
}

MPO<matrix, symm> rdm1(Model<matrix, symm> const & model, Lattice const & lattice)
{
        typedef operator_selector<matrix, symm>::type op_t;
        typedef OPTable<matrix, symm>::tag_type tag_type;
        typedef matrix::value_type value_type;

        std::shared_ptr<TagHandler<matrix, symm> > tag_handler = model.operators_table();

        const int nterm = 2;
        int pos_[nterm] = {0, 2};
        std::vector<int> pos(pos_, pos_ + nterm);

        // A
        //tag_type op1 = model.identity_matrix_tag(lattice.get_prop<symm::subcharge>("type", pos[0]));

        tag_type op1 = model.get_operator_tag("create_fill", lattice.get_prop<symm::subcharge>("type", pos[0]));
        tag_type op2 = model.get_operator_tag("destroy", lattice.get_prop<symm::subcharge>("type", pos[1]));

        //tag_type op2 = model.get_operator_tag("create_couple_up", lattice.get_prop<symm::subcharge>("type", pos[1]));
        //tag_type op3 = model.get_operator_tag("destroy_fill_couple_down", lattice.get_prop<symm::subcharge>("type", pos[2]));
        //tag_type op4 = model.get_operator_tag("destroy", lattice.get_prop<symm::subcharge>("type", pos[3]));

        tag_type ops_[nterm] = {op1, op2};
        std::vector<tag_type> ops(ops_, ops_ + nterm);

        term_descriptor<value_type> term = generate_mpo::arrange_operators(pos, ops, tag_handler);
        maquis::cout << term << std::endl;
        for (int i=0; i<term.size(); ++i)
            maquis::cout << tag_handler->get_op(term.operator_tag(i));

        std::vector<tag_type> identities, fillings;
        identities.push_back(model.identity_matrix_tag(0));
        identities.push_back(model.identity_matrix_tag(1));
        fillings.push_back(model.filling_matrix_tag(0));
        fillings.push_back(model.filling_matrix_tag(1));

        MPO<matrix, symm> mpo = generate_mpo::make_1D_mpo(pos, ops, identities, fillings, tag_handler, lattice);

        return mpo;
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

        // load state
        MPS<matrix, symm> mps;
        maquis::cout << "Loading " << parms["chkpfile"] << std::endl;
        std::string wvf = parms["chkpfile"];
        load(wvf, mps);

        MPO<matrix, symm> mpo = rdm2(model, lattice);

        maquis::traits::real_type<matrix::value_type>::type value = maquis::real(expval(mps, mpo));
        maquis::cout << "Expval is: " << value << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
