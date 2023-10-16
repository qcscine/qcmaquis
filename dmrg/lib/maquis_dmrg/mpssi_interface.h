/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef MPSSI_INTERFACE_H
#define MPSSI_INTERFACE_H

#include "maquis_dmrg.h"

namespace maquis
{
    template <class V> // real or complex
    class MPSSIInterface
    {
        public:
            // typedef for measurements
            // Intel compiler seems not to like it
            typedef maquis::meas_with_results_type<V> meas_with_results_type;
            MPSSIInterface(const std::vector<std::string>& project_names,
                           const std::vector<std::vector<int> >& states);

            ~MPSSIInterface();

            // Overlap
            V overlap(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state, bool su2u1);

            // Disabled since the interface does not need it
            // // 1-TDM
            // meas_with_results_type onetdm(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state);

            // 1-TDM, split in four spin components
            std::vector<maquis::meas_with_results_type<V> >
                onetdm_spin(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state);


            // MPS counterrotation.
            // Parameters:
            // pname: project name
            // state: state index
            // t: rotational matrix, flattened row-wise (row-major)
            // scale_inactive: inactive scaling factor
            // This function appends .rotated. to pname when saving rotated MPS
            void rotate(const std::string& pname, int state, const std::vector<V> & t, V scale_inactive, int Ms);
        private:

            // Generate 2U1 checkpoint name. rotated == true will append .rotated. to prefix
            std::string twou1_name(const std::string & pname, int state, int Ms, bool rotated = false);

            struct Impl;
            std::unique_ptr<Impl> impl_;
    };
}

#endif