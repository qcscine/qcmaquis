/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "maquis_dmrg_detail.h"
#include "dmrg/models/chem/transform_symmetry.hpp"

namespace maquis {
    namespace interface_detail {

        // Generate names for SU2U1 checkpoint files
        std::string su2u1_name(const std::string & pname, int state)
        {
            // TODO: should we check if pname and state are allowed here too with allowed_names_states()?
            // if (allowed_names_states(pname, state) == -1)
            //     throw std::runtime_error("Do not know the project name" + pname );

            std::string ret = pname + ".checkpoint_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Generate names for 2U1 checkpoint files (helper)
        std::tuple<std::string, int, int>
        twou1_name_Nup_Ndown(const std::string & pname, int state, int nel, int multiplicity, int Ms)
        {
            int Nup = (nel + Ms)/2;
            int Ndown = (nel - Ms)/2;

            int remainder = Ms;
            if (Ms == 0)
            {
                // Use 2U1 checkpoint with Ms=0 or 1 for default Ms
                remainder = nel % 2; // integer division
                Nup = nel / 2 + remainder;
                Ndown = nel / 2 - remainder;
            }

            // Use 2U1 checkpoint with Ms=0 or 1

            // int remainder = nel % 2; // integer division
            // int Nup = nel / 2 + remainder;
            // int Ndown = nel / 2 - remainder;


            std::string ret = pname + ".checkpoint_state." + std::to_string(state)
                                    + "." + std::to_string(multiplicity) + "." + std::to_string(remainder)
                                    + ".h5";
            return std::make_tuple(ret, Nup, Ndown);
        }

        std::string twou1_name(const std::string & pname, int state, int nel, int multiplicity, int Ms)
        {
            return std::get<0>(twou1_name_Nup_Ndown(pname, state, nel, multiplicity, Ms));
        }

    }

    // Set parameters required for relativistic calculations
    void prepare_relativistic(BaseParameters& parms, bool magnetic)
    {
        parms.set("symmetry", "u1dg");
        parms.set("LATTICE", "spinors");
        parms.set("CONSERVED_QUANTUMNUMBERS", "N");
        parms.set("MODEL", "relativistic_quantum_chemistry");
        parms.set("COMPLEX", "1");
        parms.set("lattice_library", "coded");
        parms.set("model_library", "coded");
        parms.set("use_compressed", "0");

        // additional relativistic options, why do we need them?
        parms.set("group_id", 8);
        parms.set("type", 0);

        // Kramer's symmetry
        if (magnetic)
            parms.set("MAGNETIC", 1);
    }

// Transforms SU2 checkpoint to 2U1 checkpoint
// Mostly copy-paste from mps_transform.cpp, but creates only one 2U1 checkpoint per state
// corresponding to the state with the highest Sz

    void transform(const std::string & pname, int state, int Ms)
    {
        // This works only for double
        // Do we really need complex? (since we don't transform double groups)
        typedef alps::numeric::matrix<double> matrix;

    #if defined(HAVE_SU2U1PG)
                typedef SU2U1PG SU2U1grp;
                typedef TwoU1PG TwoU1grp;
    #elif defined(HAVE_SU2U1)
                typedef SU2U1 SU2U1grp;
                typedef TwoU1 TwoU1grp;
    #endif

        std::string checkpoint_name = interface_detail::su2u1_name(pname, state);

        BaseParameters parms;

        if (!boost::filesystem::exists(checkpoint_name))
            throw std::runtime_error("input MPS " + checkpoint_name + " does not exist\n");

        // load source MPS
        MPS<matrix, SU2U1grp> mps;
        load(checkpoint_name, mps);

        // fetch parameters and modify symmetry
        storage::archive ar_in(checkpoint_name + "/props.h5");

        ar_in["/parameters"] >> parms;
        parms.set("init_type", "const");
    #if defined(HAVE_SU2U1PG)
        parms.set("symmetry", "2u1pg");
    #elif defined(HAVE_SU2U1)
        parms.set("symmetry", "2u1");
    #endif
        int Nup, Ndown;
        std::string twou1_checkpoint_name;
        int nel = parms["nelec"];
        int multiplicity = parms["spin"];

        // get number of up/down electrons and the checkpoint name for the 2U1 checkpoint
        std::tie(twou1_checkpoint_name, Nup, Ndown) = maquis::interface_detail::twou1_name_Nup_Ndown(pname, state, nel, multiplicity, Ms);

        parms.set("u1_total_charge1", Nup);
        parms.set("u1_total_charge2", Ndown);

        // transform MPS
        MPS<matrix, TwoU1grp> mps_out = transform_mps<matrix, SU2U1grp>()(mps, Nup, Ndown);

        save(twou1_checkpoint_name, mps_out);

        if (boost::filesystem::exists(twou1_checkpoint_name + "/props.h5"))
            boost::filesystem::remove(twou1_checkpoint_name + "/props.h5");
        boost::filesystem::copy(checkpoint_name + "/props.h5", twou1_checkpoint_name + "/props.h5");

        storage::archive ar_out(twou1_checkpoint_name + "/props.h5", "w");
        ar_out["/parameters"] << parms;
    }

}