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
            typedef typename DMRGInterface<V>::meas_with_results_type meas_with_results_type;

            MPSSIInterface(std::size_t nel,
                           const std::vector<std::size_t>& multiplicities,
                           const std::vector<std::vector<std::size_t> >& states,
                           const std::string& pname,
                           const std::vector<std::string>& mult_suffixes);

            ~MPSSIInterface();

            // Overlap
            V overlap(std::size_t bra_state, std::size_t bra_multiplicity, std::size_t ket_state, std::size_t ket_multiplicity);

            // 1-TDM
            meas_with_results_type onetdm(std::size_t bra_state, std::size_t bra_multiplicity, std::size_t ket_state, std::size_t ket_multiplicity);

        private:

            // Number of electrons. Must be same for all multiplicities (Dyson orbitals or similar not allowed yet)
            std::size_t nel_;

            // Vector with all MS2 multiplicities: 0 -- singlet, 1 -- doublet etc.
            // Can contain also same multiplicities in different entries if states have been optimised with different orbitals
            // e.g. for two state-specific singlets we have two entries in the multiplicities_ vector with 0
            std::vector<std::size_t> multiplicities_;

            //State indexes for each multiplicity
            std::vector<std::vector<std::size_t> > states_;

            // Project name
            std::string pname_;

            // Suffixes of hdf5 files for each multiplicity
            std::vector<std::string> mult_suffixes_;

            // SU2U1->2U1 transformation, takes SU2U1 checkpoint name as checkpoint_name
            void transform(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t multiplicity);

            // Generate SU2U1 checkpoint name
            std::string su2u1_name(const std::string & pname, const std::string & suffix,
                                   std::size_t state);

            // Generate 2U1 checkpoint name
            std::string twou1_name(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t multiplicity);

            // MPS counterrotation. takes 2U1 checkpoint name as a parameter!
            void rotate(const std::string& checkpoint_name);

            struct Impl;
            std::unique_ptr<Impl> impl_;
    };
}

#endif