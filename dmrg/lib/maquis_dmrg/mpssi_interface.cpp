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

#include "mpssi_interface.h"
#include "dmrg/models/chem/transform_symmetry.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/sim/matrix_types.h"

namespace maquis
{
    template <class V>
    struct MPSSIInterface<V>::Impl
    {
        // Transforms SU2 checkpoint to 2U1 checkpoint
        // Mostly copy-paste from mps_transform.cpp, but creates only one 2U1 checkpoint per state
        // corresponding to the state with the highest Sz
        void transform(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t nel, std::size_t multiplicity)
        {
#if defined(HAVE_SU2U1)
            typedef SU2U1 grp;
            typedef TwoU1 mapgrp;
#elif defined(HAVE_SU2U1PG)
            typedef SU2U1PG grp;
            typedef TwoU1PG mapgrp;
#endif
            std::string checkpoint_name = su2u1_name(pname, suffix, state);

            BaseParameters parms;

            if (!boost::filesystem::exists(checkpoint_name))
                throw std::runtime_error("input MPS " + checkpoint_name + " does not exist\n");

            // load source MPS
            MPS<matrix, grp> mps;
            load(checkpoint_name, mps);

            // fetch parameters and modify symmetry
            storage::archive ar_in(checkpoint_name + "/props.h5");

            ar_in["/parameters"] >> parms;
            parms.set("init_state", "const");
#if defined(HAVE_SU2U1)
            parms.set("symmetry", "2u1");
#elif defined(HAVE_SU2U1PG)
            parms.set("symmetry", "2u1pg");
#endif
            std::size_t Nup, Ndown;
            std::string twou1_checkpoint_name;

            // get number of up/down electrons and the checkpoint name for the 2U1 checkpoint
            std::tie(twou1_checkpoint_name, Nup, Ndown) = twou1_name_Nup_Ndown(pname, suffix, state, nel, multiplicity);

            parms.set("u1_total_charge1", Nup);
            parms.set("u1_total_charge2", Ndown);

            // transform MPS
            MPS<matrix, mapgrp> mps_out = transform_mps<matrix, grp>()(mps, Nup, Ndown);

            save(twou1_checkpoint_name, mps_out);

            if (boost::filesystem::exists(twou1_checkpoint_name + "/props.h5"))
                boost::filesystem::remove(twou1_checkpoint_name + "/props.h5");
            boost::filesystem::copy(checkpoint_name + "/props.h5", twou1_checkpoint_name + "/props.h5");

            storage::archive ar_out(twou1_checkpoint_name + "/props.h5", "w");
            ar_out["/parameters"] << parms;
        }

        std::string su2u1_name(const std::string & pname, const std::string & suffix, std::size_t state)
        {
            std::string ret = pname + "." + suffix + ".results_state." + std::to_string(state) + ".h5";
            return ret;
        }

        std::tuple<std::string, std::size_t, std::size_t>
        twou1_name_Nup_Ndown(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t nel, std::size_t multiplicity)
        {
            // Use 2U1 checkpoint with Ms=S
            std::size_t Nup = (nel + multiplicity) / 2;
            std::size_t Ndown = (nel - multiplicity) / 2;

            std::string ret = pname + "." + suffix + ".checkpoint_state." + std::to_string(state)
                                    + "." + std::to_string(multiplicity) + "." + std::to_string(Nup-Ndown)
                                    + ".h5";
            return std::make_tuple(ret, Nup, Ndown);
        }

        std::string twou1_name(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t nel, std::size_t multiplicity)
        {
            int Nup, Ndown;
            std::string ret;
            std::tie(ret, Nup, Ndown) = twou1_name_Nup_Ndown(pname, suffix, state, nel, multiplicity);
            return ret;
        }

        void rotate(const std::string & checkpoint_name)
        {

        }

        Impl() = default;
        ~Impl() = default;
    };

    template <class V>
    MPSSIInterface<V>::MPSSIInterface(std::size_t nel,
                           const std::vector<std::size_t>& multiplicities,
                           const std::vector<std::vector<std::size_t> >& states,
                           const std::string& pname,
                           const std::vector<std::string>& mult_suffixes) :
                           nel_(nel),
                           multiplicities_(multiplicities),
                           states_(states),
                           pname_(pname),
                           mult_suffixes_(mult_suffixes),
                           impl_()
    {
        assert(multiplicities_.size() == mult_suffixes_.size());
        assert(multiplicities_.size() == states_.size());

        for (std::size_t TwoS = 0; TwoS < multiplicities_.size(); TwoS++)
            for (std::size_t i; i < states[TwoS].size(); i++)
            {
                std::string chkp_name = pname_ + "." + mult_suffixes_[TwoS] + ".checkpoint_state." +
                    std::to_string(states[TwoS][i]) + ".h5";

                // Transform all checkpoints into the 2U1 representation
                transform(pname, mult_suffixes_[TwoS], states_[TwoS][i], multiplicities_[TwoS]);

                // Convert SU2 name to 2U1 name
                std::string twou1_chkp_name = twou1_name(pname_, mult_suffixes_[TwoS], states[TwoS][i], multiplicities_[TwoS]);

                // Rotate if we have more than 1 multiplicity
                if (multiplicities_.size() > 1)
                    rotate(twou1_chkp_name);
            }
    }

    template <class V>
    MPSSIInterface<V>::~MPSSIInterface() = default;

    template <class V>
    std::string MPSSIInterface<V>::su2u1_name(const std::string & pname, const std::string & suffix, std::size_t state)
    {
        return impl_->su2u1_name(pname, suffix, state);
    }

    template <class V>
    std::string MPSSIInterface<V>::twou1_name(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t multiplicity)
    {
        return impl_->twou1_name(pname, suffix, state, nel_, multiplicity);
    }

    template <class V>
    void MPSSIInterface<V>::transform(const std::string & pname, const std::string & suffix,
                                   std::size_t state, std::size_t multiplicity)
    {
        impl_->transform(pname, suffix, state, nel_, multiplicity);
    }

    template <class V>
    void MPSSIInterface<V>::rotate(const std::string & checkpoint_name) { return impl_->rotate(checkpoint_name); }

    template class MPSSIInterface<double>;
}