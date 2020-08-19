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
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/mp_tensors/mps_rotate.h"
#include "dmrg/sim/matrix_types.h"


// Internal functions

namespace maquis
{
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

    template <class V>
    struct MPSSIInterface<V>::Impl
    {

#if defined(HAVE_SU2U1PG)
        typedef SU2U1PG SU2U1grp;
        typedef TwoU1PG TwoU1grp;
#elif defined(HAVE_SU2U1)
        typedef SU2U1 SU2U1grp;
        typedef TwoU1 TwoU1grp;
#endif

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

        std::string twou1_name(std::string pname, int state, int Ms, bool rotated)
        {
            // find pname in the project names to get the correct multiplicity
            // hope this isn't too performance consuming

            int idx = allowed_names_states(pname, state);
            if (idx == -1)
                throw std::runtime_error("Do not know the project name " + pname );

            if (Ms > multiplicities_[idx])
                throw std::runtime_error("Invalid Ms in twou1_name");

            // append '.rotated' to pname if rotated == true
            if (rotated) pname += ".rotated";

            return maquis::interface_detail::twou1_name(pname, state, nel_, multiplicities_[idx], Ms);

        }

        // MPS rotation
        void rotate(const std::string & pname, int state, const std::vector<V> & t, V scale_inactive, int Ms)
        {

#if defined(HAVE_SU2U1PG)
            typedef TwoU1PG grp;
#elif defined(HAVE_SU2U1)
            typedef TwoU1 grp;
#endif
            // generate checkpoint names
            std::string checkpoint_name = twou1_name(pname, state, Ms, false);
            std::string checkpoint_name_rotated = twou1_name(pname, state, Ms, true);

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

            mps.canonize(0); // TODO: Stefan -- why is it needed?

            // check if the checkpoint has 2U1/2U1PG symmetry
            storage::archive ar_in(checkpoint_name+"/props.h5");
            BaseParameters chkp_parms;
            ar_in["/parameters"] >> chkp_parms;
            std::string sym = chkp_parms["symmetry"].str();

            if (sym.find("2u1") == std::string::npos)
                throw std::runtime_error("checkpoint for MPS rotation does not have 2U1 symmetry");

            mps_rotate::rotate_mps(mps, t_mat, scale_inactive);
            save(checkpoint_name_rotated, mps);

            // copy over props.h5 file, overwriting the old one
            if (boost::filesystem::exists(checkpoint_name_rotated + "/props.h5"))
                boost::filesystem::remove(checkpoint_name_rotated + "/props.h5");
            boost::filesystem::copy(checkpoint_name + "/props.h5", checkpoint_name_rotated + "/props.h5");

        }

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

                std::string su2u1_checkpoint_name = maquis::interface_detail::su2u1_name(pname, state);
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
                {
                    // We need to transform for all multiplicities that are needed, i.e.
                    // For each spin multiplicity, we need to have the Ms=S and Ms=min(S, otherS) for all other S
                    // so we create a list with these Ms values
                    int min_tmp = multiplicities_[i];
                    std::vector<int> mult_totransform{min_tmp};
                    for (int j = 0; j < project_names.size(); j++)
                    {
                        if (min_tmp > std::min(multiplicities_[j], min_tmp))
                        {
                            min_tmp = multiplicities_[j];
                            mult_totransform.push_back(min_tmp);
                        }
                    }

                    for (auto&& m: mult_totransform)
                        maquis::transform(pname, st, m);

                }
            }
        }
        ~Impl() = default;

        const std::vector<int> & multiplicities() const { return multiplicities_; }

        // Return multiplicity for a given pname
        int get_multiplicity(const std::string & pname)
        {
            // find the project in the project name
            auto index_itr = std::find_if(project_names_.cbegin(), project_names_.cend(), [&pname](const std::string& key)->bool{ return pname == key; });
            if (index_itr == project_names_.cend())
                throw std::runtime_error("cannot find multiplicity for project name "+pname);

            // get the corresponding multiplicity
            return multiplicities_[std::distance(project_names_.cbegin(), index_itr)];
        }

        private:

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

    };

    template <class V>
    MPSSIInterface<V>::MPSSIInterface(const std::vector<std::string>& project_names,
                           const std::vector<std::vector<int> >& states)
                           : impl_(new Impl(project_names, states))  {}
        // TODO: Implement automatic 2U1 transformation logic later


    template <class V>
    MPSSIInterface<V>::~MPSSIInterface() = default;

    template <class V>
    std::string MPSSIInterface<V>::twou1_name(const std::string & pname, int state, int Ms, bool rotated)
    {
        return impl_->twou1_name(pname, state, Ms, rotated);
    }

    // MPS rotation

    template <class V>
    void MPSSIInterface<V>::rotate(const std::string& pname, int state, const std::vector<V> & t, V scale_inactive, int Ms)
    {
        impl_->rotate(pname, state, t, scale_inactive, Ms);
    }

    // Calculate 1-TDMs
    /*
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
        parms.erase_measurements();

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
    */

    // Calculate 1-TDMs, split into two spin-components
    template <class V>
    std::vector<meas_with_results_type<V> > MPSSIInterface<V>::onetdm_spin(const std::string& bra_pname, int bra_state, const std::string& ket_pname, int ket_state)
    {

        typedef alps::numeric::matrix<V> Matrix;
        DmrgParameters parms;
        std::string ket_name, bra_name;

        // Calculate Ms according to OpenMOLCAS logic: Ms=min(S_bra, S_ket)
        int S_bra = impl_->get_multiplicity(bra_pname);
        int S_ket = impl_->get_multiplicity(ket_pname);
        int Ms = std::min(S_bra, S_ket);

        bool bra_eq_ket = (bra_pname == ket_pname) && (bra_state == ket_state);
        bool rotated = (bra_pname != ket_pname);

        ket_name = twou1_name(ket_pname, ket_state, Ms, rotated);
        bra_name = twou1_name(bra_pname, bra_state, Ms, rotated);

        storage::archive ar_in(ket_name + "/props.h5");
        ar_in["/parameters"] >> parms;


        parms.erase_measurements();

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
            std::string ket_name = maquis::interface_detail::su2u1_name(ket_pname, ket_state);
            std::string bra_name = maquis::interface_detail::su2u1_name(bra_pname, bra_state);
            return impl_->overlap_su2u1(ket_name, bra_name);
        }
        else
        {
            if ((bra_state == ket_state) && (bra_pname == ket_pname)) return (V)1.0;
            bool rotated = (bra_pname != ket_pname);

            // Actually this does not matter, but we need to provide some Ms
            const auto& multiplicities = impl_->multiplicities();
            int Ms = std::min(multiplicities[bra_state], multiplicities[ket_state]);

            std::string ket_name = twou1_name(ket_pname, ket_state, Ms, rotated);
            std::string bra_name = twou1_name(bra_pname, bra_state, Ms, rotated);
            return impl_->overlap_2u1(ket_name, bra_name);
        }
    }


    template class MPSSIInterface<double>;
}

