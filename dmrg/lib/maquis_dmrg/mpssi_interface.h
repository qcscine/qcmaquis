/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#if defined(HAVE_SU2U1PG)
            typedef SU2U1PG SU2U1grp;
            typedef TwoU1PG TwoU1grp;
#elif defined(HAVE_SU2U1)
            typedef SU2U1 SU2U1grp;
            typedef TwoU1 TwoU1grp;
#endif
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