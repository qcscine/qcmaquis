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
#include "maquis_cinterface.h"
#include <memory>
#include <string>
#include <array>
#include "maquis_dmrg.h"
#include "starting_guess.h"
#include "dmrg/utils/stdout_redirector.hpp"

std::unique_ptr<maquis::DMRGInterface<double> > interface_ptr;
DmrgParameters parms;
std::string pname;

// stdout redirector
maquis::StdoutRedirector stdout_redirect;

extern "C"
{
    typedef double V;
    void qcmaquis_interface_preinit(int nel, int L, int spin, int irrep,
                                 int* site_types, V conv_thresh, int m, int nsweeps,
                                 int* sweep_m, int nsweepm, char* project_name,
                                 bool meas_2rdm)
    {
        parms.set("symmetry", "su2u1pg");
        parms.set("nelec", nel);
        parms.set("spin", spin);
        parms.set("L", L);
        parms.set("conv_thresh", conv_thresh);
        parms.set("nsweeps", nsweeps);
        parms.set("irrep", irrep);
        parms.set("storagedir", "./tmp/");
        // for now, making things easier
        parms.set("MEASURE[1rdm]", 1);

        if (meas_2rdm) parms.set("MEASURE[2rdm]", 1);

        parms.set("max_bond_dimension", m);
        if (sweep_m != NULL)
        {
            std::string sweep_bond_dim;
            for (int i = 0; i < nsweepm; i++)
                sweep_bond_dim += std::to_string(sweep_m[i]) + ((i < nsweepm - 1) ? "," : "") ;
            parms.set("sweep_bond_dimension", sweep_bond_dim);
        }

        pname = project_name;

        // set checkpoint folder
        parms.set("chkpfile", pname + ".checkpoint_state.0.h5");
        parms.set("resultfile", pname + ".results_state.0.h5");

        if (site_types != NULL)
        {
            std::string site_types_str;
            for (int i = 0; i < L; i++)
                site_types_str += std::to_string(site_types[i]) + ((i < L - 1) ? "," : "") ;
            parms.set("site_types", site_types_str);
        }
        else
        {
            std::string site_types_str;
            for (int i = 0; i < L; i++)
                site_types_str += "0" + (i < L - 1) ? "," : "" ;
            parms.set("site_types", site_types_str);
        }

        parms.set("conv_thresh", conv_thresh);

    }

    void qcmaquis_interface_update_integrals(int* integral_indices, V* integral_values, int integral_size)
    {
        if (parms.is_set("integral_file")||parms.is_set("integrals"))
            throw std::runtime_error("updating integrals in the interface not supported yet in the FCIDUMP format");
        // set integrals
        maquis::integral_map<double> integrals;

        for (int i = 0; i < integral_size; i++)
        {
            std::array<int, 4> idx {integral_indices[4*i], integral_indices[4*i+1], integral_indices[4*i+2], integral_indices[4*i+3]};
            V value = integral_values[i];
            integrals[idx] = value;
        }

        if (!parms.is_set("integrals_binary"))
            parms.set("integrals_binary", maquis::serialize(integrals));

        // Call an integral update only if the interface has been initialised
        if (interface_ptr)
            interface_ptr->update_integrals(integrals);

    }

    void qcmaquis_interface_run_starting_guess(int nstates, char* project_name, bool do_fiedler, bool do_cideas, char* fiedler_order_string, int* hf_occupations)
    {
        // TODO: Make sure that qcmaquis_interface_preinit and _update_integrals has been called beforehand

        // if neither do_fiedler nor do_cideas are set, return

        if (!(do_fiedler || do_cideas)) return;

        std::string project_name_(project_name);

        // copy HF occupations from hf_occupations array if present
        std::vector<std::vector<int> > hf_occupations_vec;
        if (hf_occupations != nullptr)
        {
            hf_occupations_vec.reserve(nstates);
            int L = parms["L"];
            int* counter = hf_occupations;
            for (int i = 0; i < nstates; i++)
            {
                hf_occupations_vec.emplace_back(counter, counter+L);
                counter += L;
                // TODO: check for overflows
            }
        }

        maquis::StartingGuess<V> starting_guess(parms, nstates, project_name_, do_fiedler, do_cideas, hf_occupations_vec);

        if (do_fiedler)
        {
            // TODO: check length of fiedler_order_string. it must be pre-set to the correct length
            int len = strlen(fiedler_order_string);
            std::string str = starting_guess.getFiedlerOrder();
            assert(str.length() == len);
            strncpy(fiedler_order_string, str.c_str(), len);
        }

        if (do_cideas)
            starting_guess.cideas();

    }

    void qcmaquis_interface_set_nsweeps(int nsweeps)
    {
        parms.set("nsweeps", nsweeps);
    }

    void qcmaquis_interface_set_param(char* key, char* value)
    {
        parms.set(key, std::string(value));
    }

    // Start a new simulation with stored parameters
    void qcmaquis_interface_reset()
    {
        interface_ptr.reset(new maquis::DMRGInterface<double>(parms));
    }

    void qcmaquis_interface_optimize()
    {
        interface_ptr->optimize();
    }
    double qcmaquis_interface_get_energy()
    {
        return interface_ptr->energy();
    }

    void qcmaquis_interface_set_state(int state)
    {
        std::string str;
        for (int i = 0; i < state; i++)
            str += pname + ".checkpoint_state." + std::to_string(state-1) + ".h5" + ((i < state-1) ? ";" : "") ;


        parms.set("ortho_states", str);
        parms.set("n_ortho_states", state);
        parms.set("chkpfile", pname + ".checkpoint_state." + std::to_string(state) + ".h5");
        parms.set("resultfile", pname + ".results_state." + std::to_string(state) + ".h5");

        qcmaquis_interface_reset();
    }

    void qcmaquis_interface_get_1rdm(int* indices, V* values, int size)
    {
        const typename maquis::meas_with_results_type<V>& meas = interface_ptr->onerdm();
        // the size attribute is pretty much useless if we allocate the output arrays outside of the interface
        // we'll just use it to check if the size matches the size of the measurement
	//
	// When we have point group symmetry, OpenMOLCAS will allocate an array larger than our number of measurements
	// So as long as our measurements will fit into the allocated array, we're good
	// I.e. the sizes of an array do not have to match exactly
        assert(size >= meas.first.size());
        assert(size >= meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[2*i] = meas.first[i][0];
            indices[2*i+1] = meas.first[i][1];
        }
    }

    // hooray for copy-paste
    void qcmaquis_interface_get_2rdm(int* indices, V* values, int size)
    {
        const typename maquis::meas_with_results_type<V>& meas = interface_ptr->twordm();

        assert(size >= meas.first.size());
        assert(size >= meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[4*i] = meas.first[i][0];
            indices[4*i+1] = meas.first[i][1];
            indices[4*i+2] = meas.first[i][2];
            indices[4*i+3] = meas.first[i][3];

        }
    }

    // hooray for copy-paste
    void qcmaquis_interface_get_4rdm(int* indices, V* values, int size)
    {
        const typename maquis::meas_with_results_type<V>& meas = interface_ptr->fourrdm();

        assert(size >= meas.first.size());
        assert(size >= meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[8*i] = meas.first[i][0];
            indices[8*i+1] = meas.first[i][1];
            indices[8*i+2] = meas.first[i][2];
            indices[8*i+3] = meas.first[i][3];
            indices[8*i+4] = meas.first[i][4];
            indices[8*i+5] = meas.first[i][5];
            indices[8*i+6] = meas.first[i][6];
            indices[8*i+7] = meas.first[i][7];

        }
    }

    void qcmaquis_interface_measure_and_save_4rdm(int state)
    {
        std::string res_name = chem::detail::su2u1_result_name(pname, state);

        // make sure the checkpoint actually exists, otherwise we'd be measuring garbage
        // (dmrg_meas doesn't print an error yet if you try to load a non-existing checkpoint)
        if (!boost::filesystem::exists(res_name))
            throw std::runtime_error("QCMaquis checkpoint " + res_name + " does not exist. Did you optimise the wavefunction for this state?");

        // We need to remove all other measurements, so let's make a backup of current parameters and restore it when we're done
        BaseParameters parms_backup(parms);
        parms.erase_measurements();

        // reset the parameters for state state and load the checkpoint
        qcmaquis_interface_set_state(state);
        interface_ptr->measure_and_save_4rdm();

        parms = parms_backup;

    }

    void qcmaquis_interface_get_iteration_results(int* nsweeps, std::size_t* m, V* truncated_weight, V* truncated_fraction, V* smallest_ev)
    {
        // Get iteration results from the last sweep
        results_collector& iter = interface_ptr->get_iteration_results();

        // iter contains results, one element per microiteration
        const std::vector<boost::any>& m_vec = iter["BondDimension"].get();
        const std::vector<boost::any>& tw_vec = iter["TruncatedWeight"].get();
        const std::vector<boost::any>& tf_vec = iter["TruncatedFraction"].get();
        const std::vector<boost::any>& ev_vec = iter["SmallestEV"].get();

        // We return the sum of these values for the last sweep
        // this should be done with transform_reduce

        *m = 0; for (auto&& m_ : m_vec) *m += boost::any_cast<std::size_t>(m_);
        *truncated_weight = 0; for (auto&& tw_ : tw_vec) *truncated_weight += boost::any_cast<V>(tw_);
        *truncated_fraction = 0; for (auto&& tf_ : tf_vec) *truncated_fraction += boost::any_cast<V>(tf_);
        *smallest_ev = 0; for (auto&& ev_ : ev_vec) *smallest_ev += boost::any_cast<V>(ev_);
        *nsweeps = interface_ptr->get_last_sweep();
    }

    double qcmaquis_interface_get_overlap(char* filename)
    {
        return interface_ptr->overlap(filename);
    }

    void qcmaquis_interface_stdout(char* filename)
    {
        stdout_redirect.set_filename(filename);
    }

    void qcmaquis_interface_restore_stdout()
    {
        stdout_redirect.restore();
    }
}
