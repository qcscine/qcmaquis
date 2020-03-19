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
        void transform(const std::string & pname, int state, bool Ms_equal_s=false)
        {
#if defined(HAVE_SU2U1PG)
            typedef SU2U1PG grp;
            typedef TwoU1PG mapgrp;
#elif defined(HAVE_SU2U1)
            typedef SU2U1 grp;
            typedef TwoU1 mapgrp;
#endif
            std::string checkpoint_name = su2u1_name(pname, state);

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
            std::tie(twou1_checkpoint_name, Nup, Ndown) = twou1_name_Nup_Ndown(pname, state, nel, multiplicity, Ms_equal_s);

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
        std::string su2u1_name(const std::string & pname, int state)
        {
            // TODO: should we check if pname and state are allowed here too with allowed_names_states()?
            // if (allowed_names_states(pname, state) == -1)
            //     throw std::runtime_error("Do not know the project name" + pname );

            std::string ret = pname + ".checkpoint_state." + std::to_string(state) + ".h5";
            return ret;
        }

        // Generate names for 2U1 checkpoint files
        std::tuple<std::string, int, int>
        twou1_name_Nup_Ndown(const std::string & pname, int state, int nel, int multiplicity, bool Ms_equal_s = false)
        {
            int Nup, Ndown;
            std::string ret;

            if (Ms_equal_s)
            {
                // Use 2U1 checkpoint with Ms=S
                Nup = (nel + multiplicity) / 2;
                Ndown = (nel - multiplicity) / 2;
                ret = pname + ".checkpoint_state." + std::to_string(state)
                            + "." + std::to_string(multiplicity) + "." + std::to_string(Nup-Ndown)
                            + ".h5";
            }
            else
            {
                // Use 2U1 checkpoint with Ms=0 or 1
                int remainder = nel % 2; // integer division
                Nup = nel / 2 + remainder;
                Ndown = nel / 2 - remainder;
                ret = pname + ".checkpoint_state." + std::to_string(state)
                            + "." + std::to_string(multiplicity) + "." + std::to_string(remainder)
                            + ".h5";
            }


            return std::make_tuple(ret, Nup, Ndown);
        }

        std::string twou1_name(const std::string & pname, int state, bool Ms_equal_s = false)
        {
            // find pname in the project names to get the correct multiplicity
            // hope this isn't too performance consuming

            int idx = allowed_names_states(pname, state);
            if (idx == -1)
                throw std::runtime_error("Do not know the project name " + pname );

            int Nup, Ndown;
            std::string ret;
            std::tie(ret, Nup, Ndown) = twou1_name_Nup_Ndown(pname, state, nel_, multiplicities_[idx], Ms_equal_s);
            return ret;
        }

        // MPS rotation
        void rotate(const std::string & checkpoint_name, const std::vector<V> & t, V scale_inactive)
        {

#if defined(HAVE_SU2U1PG)
            typedef TwoU1PG grp;
#elif defined(HAVE_SU2U1)
            typedef TwoU1 grp;
#endif
            // convert t to alps::matrix

            // check if the length of the vector is a square number
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

            if (sym.find("2u1") == std::string::npos)
                throw std::runtime_error("checkpoint for MPS rotation does not have 2U1 symmetry");

            mps_rotate::rotate_mps(mps, t_mat, scale_inactive);
            save(checkpoint_name, mps);
        }

        // Check if the project name and the state index are allowed
        // (i.e. the project prefix and the state is present in the vectors with which the class has been initialised)
        // If project name or state is not found, returns -1
        // otherwise returns the index of the project name (we don't care about the state index for now)
        // --- hope searching a vector of strings isn't too performance consuming ---
        int allowed_names_states(const std::string & pname, int state)
        {
            // check if pname is allowed
            auto index_itr = std::find_if(project_names_.cbegin(), project_names_.cend(), [&pname](const std::string& key)->bool{ return pname == key; });
            if (index_itr == project_names_.cend())
                return -1;

            int idx = std::distance(project_names_.cbegin(), index_itr);

            // check if the state is found
            auto index_st_itr = std::find_if(states_[idx].cbegin(), states_[idx].cend(), [&state](const int key)->bool{ return state == key; });
            if (index_st_itr == states_[idx].cend())
                return -1;

            return idx;
        }

        // State indexes for each project. E.g. if states_[0] is {0,3}, we are interested in states 0 and 3 of the first project
        // Each project has a different name, a typical use would be to use one project name for each multiplicity
        // or more general, results of a single DMRGSCF calculation
        std::vector<std::vector<int> > states_;

        // Suffixes of hdf5 files for each multiplicity
        std::vector<std::string> project_names_;

        // All multiplicities
        std::vector<int> multiplicities_;

        // Number of electrons, for now only one number of electrons supported for all projects. (i.e. no Dyson orbitals)
        int nel_;

        // Constructor that takes project names and states to get the multiplicities
        Impl(const std::vector<std::string> & project_names, const std::vector<std::vector<int> >& states)
            : project_names_(project_names), nel_(-1), multiplicities_(states.size()), states_(states)
        {
            assert(project_names.size() == states.size());

            // Find out about the spin multiplicity required in the 2U1 filename
            // by loading the parameters from the SU2U1 file

            for (int i = 0; i < project_names.size(); i++)
            {
                // Load the first state of each project group and get its number of electrons and spin.
                // We will not check if the subsequent states have the same spin for now, but if you want to implement if for better
                // error safety, feel free.
                const std::string& pname = project_names[i];
                assert(states[i].size() > 0); // make sure we do not have empty state containers
                int state = states[i][0];

                std::string su2u1_checkpoint_name = su2u1_name(pname, state);
                BaseParameters parms;
                if (!boost::filesystem::exists(su2u1_checkpoint_name))
                    throw std::runtime_error("SU2U1 MPS checkpoint " + su2u1_checkpoint_name + " is required but does not exist\n");

                storage::archive ar_in(su2u1_checkpoint_name + "/props.h5");

                ar_in["/parameters"] >> parms; // TODO: check if /parameters/nelec and /parameters/spin can be read directly and if it saves time

                int nel = parms["nelec"];
                if ((nel_ != -1) && (nel_ != nel))
                    throw std::runtime_error("Different electron numbers in different project groups: This is not supported in MPSSI yet.");
                nel_ = nel;
                multiplicities_[i] = parms["spin"];

                // transform all checkpoints to 2U1 point group

                for (auto&& st: states[i])
                    transform(pname, st);
            }
        }
        ~Impl() = default;

    };

    template <class V>
    MPSSIInterface<V>::MPSSIInterface(const std::vector<std::string>& project_names,
                           const std::vector<std::vector<int> >& states)
                           : impl_(new Impl(project_names, states))  {}
        // TODO: Implement automatic 2U1 transformation logic later


    template <class V>
    MPSSIInterface<V>::~MPSSIInterface() = default;

    template <class V>
    std::string MPSSIInterface<V>::su2u1_name(const std::string & pname, int state)
    {
        return impl_->su2u1_name(pname, state);
    }

    template <class V>
    std::string MPSSIInterface<V>::twou1_name(const std::string & pname, int state, bool Ms_equal_s)
    {
        return impl_->twou1_name(pname, state, Ms_equal_s);
    }

    // SU2U1->2U1 transformation
    template <class V>
    void MPSSIInterface<V>::transform(const std::string & pname, int state, bool Ms_equal_s)
    {
        impl_->transform(pname, state, Ms_equal_s);
    }

    // MPS rotation
    template <class V>
    void MPSSIInterface<V>::rotate(const std::string & checkpoint_name, const std::vector<V> & t, V scale_inactive)
    {
        impl_->rotate(checkpoint_name, t, scale_inactive);
    }

    template <class V>
    void MPSSIInterface<V>::rotate(const std::string& pname, int state, const std::vector<V> & t, V scale_inactive)
    {
        impl_->rotate(twou1_name(pname, state), t, scale_inactive);
    }

    // Calculate 1-TDMs
    template <class V>
    meas_with_results_type<V> MPSSIInterface<V>::onetdm(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state)
    {
        typedef alps::numeric::matrix<V> Matrix;

        DmrgParameters parms;

        std::string ket_name, bra_name;

        bool bra_eq_ket = (bra_pname == ket_pname) && (bra_state == ket_state);
        if (bra_pname == ket_pname) // if we have the same multiplicity, use SU2U1
        {
            ket_name = su2u1_name(ket_pname, ket_state);
            bra_name = su2u1_name(bra_pname, bra_state);
#if defined(HAVE_SU2U1PG)
            parms.set("symmetry", "su2u1pg");
#elif defined(HAVE_SU2U1)
            parms.set("symmetry", "su2u1");
#endif
        }
        else // otherwise, 2U1
        {
            ket_name = twou1_name(ket_pname, ket_state);
            bra_name = twou1_name(bra_pname, bra_state);
#if defined(HAVE_SU2U1PG)
            parms.set("symmetry", "2u1pg");
#elif defined(HAVE_SU2U1)
            parms.set("symmetry", "2u1");
#endif
        }

        // load all parameters from HDF5, but remove measurements
        storage::archive ar_in(ket_name + "/props.h5");
        ar_in["/parameters"] >> parms;
        parms.erase_substring("MEASURE");

        parms.set("chkpfile", ket_name);
        if (bra_eq_ket) // run 1-RDM measurement if bra == ket
            parms.set("MEASURE[1rdm]", "1");
        else
            parms.set("MEASURE[trans1rdm]", bra_name);

        // run measurement
        maquis::DMRGInterface<double> interface(parms);
        interface.measure();
        if (bra_eq_ket)
            return interface.measurements().at("oneptdm");
        return interface.measurements().at("transition_oneptdm");
    }

    // Calculate 1-TDMs, split into two spin-components
    template <class V>
    std::vector<meas_with_results_type<V> > MPSSIInterface<V>::onetdm_spin(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state)
    {

        typedef alps::numeric::matrix<V> Matrix;
        DmrgParameters parms;
        std::string ket_name, bra_name;

        // For some reason for bra==ket MPSSI is using Ms==0 or 1 instead of Ms == S

        bool bra_eq_ket = (bra_pname == ket_pname) && (bra_state == ket_state);

        ket_name = twou1_name(ket_pname, ket_state, bra_eq_ket);
        bra_name = twou1_name(bra_pname, bra_state, bra_eq_ket);

        if (bra_eq_ket)
            transform(bra_pname, bra_state, bra_eq_ket);

        storage::archive ar_in(ket_name + "/props.h5");
        ar_in["/parameters"] >> parms;


        parms.erase_substring("MEASURE");

#if defined(HAVE_SU2U1PG)
        parms.set("symmetry", "2u1pg");
#elif defined(HAVE_SU2U1)
        parms.set("symmetry", "2u1");
#endif

        parms.set("chkpfile", ket_name);
        if (bra_eq_ket) // run 1-RDM measurement if bra == ket
        {
            parms.set("MEASURE[1rdm_aa]", "1");
            // ab and ba are not needed in MOLCAS
            // parms.set("MEASURE[1rdm_ab]", "1");
            // parms.set("MEASURE[1rdm_ba]", "1");
            parms.set("MEASURE[1rdm_bb]", "1");
        }
        else
        {
            parms.set("MEASURE[trans1rdm_aa]", bra_name);
            // parms.set("MEASURE[trans1rdm_ab]", bra_name);
            // parms.set("MEASURE[trans1rdm_ba]", bra_name);
            parms.set("MEASURE[trans1rdm_bb]", bra_name);
        }

        // run measurement
        maquis::DMRGInterface<double> interface(parms);
        interface.measure();

        std::vector<meas_with_results_type> ret;
        ret.reserve(2);

        const results_map_type<V>& meas = interface.measurements();

        if (bra_eq_ket)
        {
            ret.push_back(meas.at("oneptdm_aa"));
            // ret.push_back(meas.at("oneptdm_ab"));
            // ret.push_back(meas.at("oneptdm_ba"));
            ret.push_back(meas.at("oneptdm_bb"));

        }
        else
        {
            ret.push_back(meas.at("transition_oneptdm_aa"));
            // ret.push_back(meas.at("transition_oneptdm_ab"));
            // ret.push_back(meas.at("transition_oneptdm_ba"));
            ret.push_back(meas.at("transition_oneptdm_bb"));
        }

        return ret;
    }


    template <class V>
    V MPSSIInterface<V>::overlap(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state, bool su2u1)
    {
        if ((bra_pname == ket_pname) && su2u1)
        {
            if (bra_state == ket_state) return (V)1.0;
            std::string ket_name = su2u1_name(ket_pname, ket_state);
            std::string bra_name = su2u1_name(bra_pname, bra_state);
            return impl_->overlap_su2u1(ket_name, bra_name);
        }
        else
        {
            if ((bra_state == ket_state) && (bra_pname == ket_pname)) return (V)1.0;
            std::string ket_name = twou1_name(ket_pname, ket_state);
            std::string bra_name = twou1_name(bra_pname, bra_state);
            return impl_->overlap_2u1(ket_name, bra_name);
        }
    }


    template class MPSSIInterface<double>;
}

