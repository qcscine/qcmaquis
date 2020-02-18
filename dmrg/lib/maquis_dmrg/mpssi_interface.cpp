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
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "dmrg/sim/matrix_types.h"

// Internal functions
namespace detail
{
    // checks if x is an integer square
    template <class I>
    bool is_square(I x)
    {
        // we only need this for integer types, but floating point should also work ok
        // if we use floating point, then the return of the sqrt should be of the same type
        typedef typename std::conditional<std::is_floating_point<I>::value, I, double>::type T;

        T sr = sqrt(x);

        // If square root is an integer
        return ((sr - floor(sr)) == 0);
    }
}

namespace maquis
{
    template <class V>
    struct MPSSIInterface<V>::Impl
    {
        typedef alps::numeric::matrix<V> Matrix;
        // Overlap calculation
        V overlap_2u1(const std::string& bra_checkpoint, const std::string& ket_checkpoint)
        {
            MPS<Matrix, TwoU1grp> bra, ket;
            load(bra_checkpoint, bra);
            load(ket_checkpoint, ket);
            return (V) ::overlap(bra, ket);
        }

        V overlap_su2u1(const std::string& bra_checkpoint, const std::string& ket_checkpoint)
        {
            MPS<matrix, SU2U1grp> bra, ket;
            load(bra_checkpoint, bra);
            load(ket_checkpoint, ket);
            return (V) ::overlap(bra, ket);
        }

        // Transforms SU2 checkpoint to 2U1 checkpoint
        // Mostly copy-paste from mps_transform.cpp, but creates only one 2U1 checkpoint per state
        // corresponding to the state with the highest Sz
        void transform(const std::string & pname, const std::string & suffix,
                                   int state, int nel, int multiplicity)
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
            MPS<matrix, SU2U1grp> mps;
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
            int Nup, Ndown;
            std::string twou1_checkpoint_name;

            // get number of up/down electrons and the checkpoint name for the 2U1 checkpoint
            std::tie(twou1_checkpoint_name, Nup, Ndown) = twou1_name_Nup_Ndown(pname, suffix, state, nel, multiplicity);

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

        // Generate names for SU2U1 checkpoint files
        std::string su2u1_name(const std::string & pname, const std::string & suffix, int state)
        {
            std::string ret = pname + "." + suffix + ".results_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Generate names for 2U1 checkpoint files
        std::tuple<std::string, int, int>
        twou1_name_Nup_Ndown(const std::string & pname, const std::string & suffix,
                                   int state, int nel, int multiplicity)
        {
            // Use 2U1 checkpoint with Ms=S
            int Nup = (nel + multiplicity) / 2;
            int Ndown = (nel - multiplicity) / 2;

            std::string ret = pname + "." + suffix + ".checkpoint_state." + std::to_string(state)
                                    + "." + std::to_string(multiplicity) + "." + std::to_string(Nup-Ndown)
                                    + ".h5";
            return std::make_tuple(ret, Nup, Ndown);
        }

        std::string twou1_name(const std::string & pname, const std::string & suffix,
                                   int state, int nel, int multiplicity)
        {
            int Nup, Ndown;
            std::string ret;
            std::tie(ret, Nup, Ndown) = twou1_name_Nup_Ndown(pname, suffix, state, nel, multiplicity);
            return ret;
        }

        // MPS rotation
        void rotate(const std::string & checkpoint_name, const std::vector<V> & t, V scale_inactive)
        {

#if defined(HAVE_SU2U1)
            typedef TwoU1 grp;
#elif defined(HAVE_SU2U1PG)
            typedef TwoU1PG grp;
#endif
            // convert t to alps::matrix

            // check if the length of the vector is a triangular number,
            // a necessary condition for the size of a flattened full lower-triangular matrix
            assert(detail::is_square(t.size()));

            // extract matrix dimension from the vector size
            int dim = floor(sqrt(t.size()));

            Matrix t_mat(dim,dim);

            int idx = 0;
            for (int i = 0; i < dim; i++)
                for (int j = 0; j < dim; j++)
                    t_mat(i,j) = t[idx++];


            MPS<Matrix, grp> mps;
            load(checkpoint_name, mps);

            // check if the checkpoint has 2U1/2U1PG symmetry
            storage::archive ar_in(checkpoint_name+"/props.h5");
            BaseParameters chkp_parms;
            ar_in["/parameters"] >> chkp_parms;
            std::string sym = chkp_parms["symmetry"].str();

            if (sym.find("2u1") != std::string::npos)
                throw std::runtime_error("checkpoint for MPS rotation does not have 2U1 symmetry");

            mps_rotate::rotate_mps(mps, t_mat, scale_inactive);
            save(checkpoint_name, mps);
        }

        Impl() = default;
        ~Impl() = default;
    };

    template <class V>
    MPSSIInterface<V>::MPSSIInterface(int nel,
                           const std::vector<int>& multiplicities,
                           const std::vector<std::vector<int> >& states,
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

        // TODO: Implement automatic 2U1 transformation logic later
        // for (int TwoS = 0; TwoS < multiplicities_.size(); TwoS++)
        //     for (int i; i < states[TwoS].size(); i++)
        //     {
                // std::string chkp_name = pname_ + "." + mult_suffixes_[TwoS] + ".checkpoint_state." +
                //     std::to_string(states[TwoS][i]) + ".h5";

                // Transform all checkpoints into the 2U1 representation
                // transform(pname, mult_suffixes_[TwoS], states_[TwoS][i], multiplicities_[TwoS]);

                // // Convert SU2 name to 2U1 name
                // std::string twou1_chkp_name = twou1_name(pname_, mult_suffixes_[TwoS], states[TwoS][i], multiplicities_[TwoS]);
            // }
    }

    template <class V>
    MPSSIInterface<V>::~MPSSIInterface() = default;

    template <class V>
    std::string MPSSIInterface<V>::su2u1_name(const std::string & pname, const std::string & suffix, int state)
    {
        return impl_->su2u1_name(pname, suffix, state);
    }

    template <class V>
    std::string MPSSIInterface<V>::twou1_name(const std::string & pname, const std::string & suffix,
                                   int state, int multiplicity)
    {
        return impl_->twou1_name(pname, suffix, state, nel_, multiplicity);
    }

    // SU2U1->2U1 transformation
    template <class V>
    void MPSSIInterface<V>::transform(const std::string & pname, const std::string & suffix,
                                   int state, int multiplicity)
    {
        impl_->transform(pname, suffix, state, nel_, multiplicity);
    }

    // MPS rotation
    template <class V>
    void MPSSIInterface<V>::rotate(const std::string & checkpoint_name, const std::vector<V> & t, V scale_inactive)
    {
        return impl_->rotate(checkpoint_name, t, scale_inactive);
    }

    // Calculate 1-TDMs
    template <class V>
    meas_with_results_type<V> MPSSIInterface<V>::onetdm(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity)
    {
        typedef alps::numeric::matrix<V> Matrix;

        DmrgParameters parms;
        const chem::integral_map<typename Matrix::value_type> fake_integrals = { { { 1, 1, 1, 1 },   0.0 } };

        parms.set("integrals_binary", chem::serialize(fake_integrals));

        std::string ket_name, bra_name;
        if (bra_multiplicity == ket_multiplicity) // if we have the same multiplicity, use SU2U1
        {
            ket_name = su2u1_name(pname_, mult_suffixes_[ket_multiplicity], ket_state);
            bra_name = su2u1_name(pname_, mult_suffixes_[bra_multiplicity], bra_state);
#if defined(HAVE_SU2U1)
            parms.set("symmetry", "su2u1");
#elif defined(HAVE_SU2U1PG)
            parms.set("symmetry", "su2u1pg");
#endif
        }
        else // otherwise, 2U1
        {
            ket_name = twou1_name(pname_, mult_suffixes_[ket_multiplicity], ket_state, ket_multiplicity);
            bra_name = twou1_name(pname_, mult_suffixes_[bra_multiplicity], bra_state, bra_multiplicity);
#if defined(HAVE_SU2U1)
            parms.set("symmetry", "2u1");
#elif defined(HAVE_SU2U1PG)
            parms.set("symmetry", "2u1pg");
#endif
        }
        parms.set("chkpfile", ket_name);
        if (bra_state == ket_state) // run 1-RDM measurement if bra == ket
            parms.set("MEASURE[1rdm]", "1");
        else
            parms.set("MEASURE[trans1rdm]", bra_name);

        // run measurement
        maquis::DMRGInterface<double> interface(parms);
        interface.measure();
        if (bra_state == ket_state)
            return interface.measurements().at("oneptdm");
        return interface.measurements().at("transition_oneptdm");
    }

    template <class V>
    V MPSSIInterface<V>::overlap(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity, bool su2u1)
    {
        if ((bra_multiplicity == ket_multiplicity) && su2u1)
        {
            if (bra_state == ket_state) return (V)1.0;
            std::string ket_name = su2u1_name(pname_, mult_suffixes_[ket_multiplicity], ket_state);
            std::string bra_name = su2u1_name(pname_, mult_suffixes_[bra_multiplicity], bra_state);
            return impl_->overlap_su2u1(ket_name, bra_name);
        }
        else
        {
            if ((bra_state == ket_state) && (bra_multiplicity == ket_multiplicity)) return (V)1.0;
            std::string ket_name = twou1_name(pname_, mult_suffixes_[ket_multiplicity], ket_state, ket_multiplicity);
            std::string bra_name = twou1_name(pname_, mult_suffixes_[bra_multiplicity], bra_state, bra_multiplicity);
            return impl_->overlap_2u1(ket_name, bra_name);
        }
    }


    template class MPSSIInterface<double>;
}