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
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/sim/matrix_types.h"
#include "dmrg/block_matrix/indexing.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/mp_tensors/mps_initializers.h"

#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps_join.h"

#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/chem/util.h"

#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#elif defined(USE_U1DG)
typedef U1DG grp;
#endif

int main(int argc, char ** argv)
{
    try {
        //if (argc != 2) {
        //    std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
        //    return 1;
        //}

        int L = 4;
        int Nup = 2;
        int Ndown = 2;
        grp::subcharge irrep = 0;

        BaseParameters parms = chem::detail::set_2u1_parameters(L, Nup, Ndown);
        parms.set("init_bond_dimension", 1000);
        parms.set("site_types", "0,0,0,0");

        default_mps_init<matrix, grp> mpsinit(parms,
                                            chem::detail::make_2u1_site_basis<matrix, grp>(L, Nup, Ndown, parms["site_types"]),
                                            chem::detail::make_2u1_initc<grp>(Nup, Ndown, irrep), parms["site_types"]);
        MPS<matrix, grp> mps(L, mpsinit);
        save("mps.h5", mps);

        MPS<matrix, grp> mps1(L, mpsinit);
        save("mps1.h5", mps1);
        MPS<matrix, grp> mps2(L, mpsinit);
        save("mps2.h5", mps2);

        Lattice lat(parms);
        Model<matrix, grp> model(lat, parms);

        typedef Lattice::pos_t pos_t;
        typedef OPTable<matrix, grp>::tag_type tag_type;

        std::vector<pos_t> positions;
        std::vector<tag_type> operators;
        std::vector<tag_type> ident, fill;

        positions.push_back(1);
        positions.push_back(2);
        operators.push_back(model.get_operator_tag("destroy_up", 0));
        operators.push_back(model.get_operator_tag("create_up", 0));
        ident.push_back(model.identity_matrix_tag(0));
        fill.push_back(model.filling_matrix_tag(0));

        MPO<matrix, grp> mpo = generate_mpo::make_1D_mpo(positions, operators, ident, fill, model.operators_table(), lat);

        grp::charge delta = grp::IdentityCharge;
        MPS<matrix, grp> pmps(lat.size());
        for (int p = 0; p < lat.size(); ++p) {
            maquis::cout << "\n****************** site " << p << " ****************\n";
            pmps[p] =  mpo_times_mps(mpo[p], mps[p], delta);
        }
        save("pmps.h5", pmps);

        MPS<matrix, grp> amps = join(pmps, mps);
        save("amps.h5", amps);
        maquis::cout << mps[1] << std::endl;
        maquis::cout << pmps[1] << std::endl;
        maquis::cout << amps[1] << std::endl;


    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
