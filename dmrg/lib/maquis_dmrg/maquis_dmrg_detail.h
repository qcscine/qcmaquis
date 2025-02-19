/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef MAQUIS_DMRG_DETAIL_H
#define MAQUIS_DMRG_DETAIL_H

#include "integral_interface.h"
#include "dmrg/utils/results_collector.h"
#include "dmrg/utils/DmrgParameters.h"

namespace maquis {
    using chem::integral_map;
    using chem::serialize;
    // Constants for integral_map template specialisation, whether one should use relativistic or nonrelativistic integrals
    namespace integrals {
        const bool relativistic = true;
        const bool nonrelativistic = false;
    }

    namespace interface_detail {
        // Generate names for SU2U1 checkpoint files
        std::string su2u1_name(const std::string & pname, int state);

        // Generate names for 2U1 checkpoint files
        // Generates only the checkpoint with minimum Ms, as consistently used in the interface
        // returns a tuple containing the filename, Nup and Ndown electrons
        // If Ms is set to 0, the Ms closest to 0 (i.e. 0 or 1) will be chosen
        // (to allow for default values where we do not care for Ms)
        std::tuple<std::string, int, int> twou1_name_Nup_Ndown(const std::string & pname, int state, int nel, int multiplicity, int Ms=0);

        // same as above, but returns only the filename
        std::string twou1_name(const std::string & pname, int state, int nel, int multiplicity, int Ms=0);

        // SU2U1 result name file
        // TODO: Similar functions for (SU2U1 and 2U1) checkpoint names are found in mpssi_interface.cpp, maybe they should be moved here?
        inline std::string su2u1_result_name(const std::string& pname, int state)
        {
            std::string ret = pname + ".results_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Result file name for transition 3-RDM measurements
        inline std::string trans3rdm_result_name(const std::string& pname, int state, int bra_state)
        {
            std::string ret = pname + ".trans3rdm." + std::to_string(state) + "_" + std::to_string(bra_state) + ".h5";
            return ret;
        }

        // Convert project name to working directory by striping last characters
        inline std::string pname2workdir(const std::string& pname) {
          // Find the last occurrence of '/'
          size_t lastSlashPos = pname.find_last_of('/');

          // If '/' is found, return the substring up to and including the last '/'
          if (lastSlashPos != std::string::npos) {
            return pname.substr(0, lastSlashPos + 1);
          }

          // If no '/' is found, return the original string (or handle accordingly)
          return pname;
        }

        // Extacts state number from chkpfile string
        // Expects {system}.{state_num}.h5
        int get_state_number(const std::string& chkpfile) {
            // Find positions of the last two periods
          size_t lastDot = chkpfile.find_last_of('.'); // Last '.'
          size_t secondLastDot = chkpfile.find_last_of('.', lastDot - 1); // Second last '.'
          
          int state_num = 0;

        if (secondLastDot != std::string::npos && lastDot != std::string::npos) {
            std::string numberStr = chkpfile.substr(secondLastDot + 1, lastDot - secondLastDot - 1); // Extract substring
            int stateNumber = std::stoi(numberStr); // Convert to integer
            std::cout << "State number: " << stateNumber << std::endl;
            
        } else {
            std::cout << "State number could not be extracted, using 0" << std::endl;
        }
        return state_num;
        }

    // Set parameters required for relativistic calculation
    void prepare_relativistic(BaseParameters& parms, bool magnetic = false);

    // Transforms SU2 checkpoint to 2U1 checkpoint
    void transform(const std::string & pname, int state, int Ms=0);

}


#endif
