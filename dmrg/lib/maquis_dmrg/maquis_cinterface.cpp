/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#include "maquis_cinterface.h"
#include <memory>
#include <string>
#include <array>
#include <regex>
#include "maquis_dmrg.h"
#include "starting_guess.h"
#include "dmrg/utils/stdout_redirector.hpp"
#include "dmrg/models/measurements/measurements_details.h" // for 4-RDM functions
                                                           
// For CASPT2
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"

std::unique_ptr<maquis::DMRGInterface<double> > interface_ptr;
DmrgParameters parms;
std::string pname;

// stdout redirector
maquis::StdoutRedirector stdout_redirect;

extern "C"
{
    typedef double V;
    void qcmaquis_interface_preinit(int nel, int L, int spin, int irrep,
                                 const int* site_types, V conv_thresh, int m, int nsweeps,
                                 const int* sweep_m, int nsweepm, const char* project_name,
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

    void qcmaquis_interface_preinit_checkpoint(const char* checkpoint_name)
    {
        // extract project name from checkpoint
        std::regex r("^(.+)\\.checkpoint");
        std::smatch m;
        std::string name_str(checkpoint_name);
        std::regex_search(name_str, m, r);
        // match found
        if (m.size() <= 1)
            throw std::runtime_error("Cannot deduce project name from the checkpoint name");

        pname = m[1];
        std::string props_name = std::string(checkpoint_name) + "/props.h5";
        if (!boost::filesystem::exists(props_name))
            throw std::runtime_error("Filename " + props_name + " cannot be found.");

        storage::archive props(props_name);
        props["/parameters"] >> parms;

    }

    void qcmaquis_interface_update_integrals(const int* integral_indices, const V* integral_values, int integral_size)
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

    void qcmaquis_interface_run_starting_guess(int nstates, const char* project_name, bool do_fiedler, bool do_cideas, char* fiedler_order_string, int* hf_occupations)
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

    void qcmaquis_interface_set_param(const char* key, const char* value)
    {
        parms.set(key, std::string(value));
    }

    void qcmaquis_interface_remove_param(const char* key)
    {
        parms.erase(key);
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
            str += pname + ".checkpoint_state." + std::to_string(i) + ".h5" + ((i < state-1) ? " " : "") ;


        parms.set("ortho_states", str);
        parms.set("n_ortho_states", state);
        parms.set("chkpfile", maquis::interface_detail::su2u1_name(pname, state));
        parms.set("resultfile", maquis::interface_detail::su2u1_result_name(pname, state));

        qcmaquis_interface_reset();
    }

    // Generic function to request either 1-RDM or spin-DM
    // to avoid copy-paste between get_1rdm and get_spdm
    // meas == interface_ptr->onerdm(): 1-RDM
    // meas == interface_ptr->onespdm(): spin-DM
    void qcmaquis_interface_get_generic1rdm(const typename maquis::meas_with_results_type<V>& meas, int* indices, V* values, int size)
    {
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

    void qcmaquis_interface_get_1rdm(int* indices, V* values, int size)
    {
        qcmaquis_interface_get_generic1rdm(interface_ptr->onerdm(), indices, values, size);
    }

    void qcmaquis_interface_get_spdm(int* indices, V* values, int size)
    {
        qcmaquis_interface_get_generic1rdm(interface_ptr->onespdm(), indices, values, size);
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
    void qcmaquis_interface_get_3rdm(int* indices, V* values, int size)
    {
        const typename maquis::meas_with_results_type<V>& meas = interface_ptr->threerdm();

        assert(size >= meas.first.size());
        assert(size >= meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[6*i] = meas.first[i][0];
            indices[6*i+1] = meas.first[i][1];
            indices[6*i+2] = meas.first[i][2];
            indices[6*i+3] = meas.first[i][3];
            indices[6*i+4] = meas.first[i][4];
            indices[6*i+5] = meas.first[i][5];

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

    #define measure_and_save_rdm(N, state) \
    std::string res_name = maquis::interface_detail::su2u1_result_name(pname, state); \
    BaseParameters meas_parms = parms.measurements(); \
    parms.erase_measurements(); \
    parms.set("MEASURE[" #N "rdm]", 1); \
    qcmaquis_interface_set_state(state); \
    interface_ptr->measure_and_save_ ## N ## rdm(); \
    parms.erase_measurements(); \
    parms << meas_parms

    void qcmaquis_interface_measure_and_save_4rdm(int state)
    {
        measure_and_save_rdm(4, state);
    }
    void qcmaquis_interface_measure_and_save_3rdm(int state)
    {
        measure_and_save_rdm(3, state);
    }
    #undef measure_and_save_rdm

    // this one does not use the macro above as it's a bit too complicated
    void qcmaquis_interface_measure_and_save_trans3rdm(int state, int bra_state)
    {

        // change the result file name to a special name that contains state numbers
        // because in an ordinary HDF5 result file it is only possible to have one transition 3-RDM
        // but we need one for each bra_state
        std::string bra_chkp = maquis::interface_detail::su2u1_name(pname, bra_state);

        std::string rfile = maquis::interface_detail::trans3rdm_result_name(pname, state, bra_state);
        std::string old_rfile;
        if (parms.is_set("resultfile"))
            old_rfile = parms["resultfile"].str();

        if (!boost::filesystem::exists(bra_chkp))
            throw std::runtime_error("QCMaquis checkpoint " + bra_chkp +
                " does not exist. Did you optimise the wavefunction for this state?");
        BaseParameters meas_parms = parms.measurements();
        parms.erase_measurements();
        parms.set("MEASURE[trans3rdm]", bra_chkp);
        qcmaquis_interface_set_state(state);
        parms["resultfile"] = rfile;

        interface_ptr->measure_and_save_trans3rdm(bra_chkp);
        parms.erase_measurements();
        parms << meas_parms;

        // restore old result file name
        if (old_rfile.empty())
            parms.erase("resultfile");
        else
            parms["resultfile"] = old_rfile;
    }

    void qcmaquis_interface_get_iteration_results(int* nsweeps, std::size_t* m, V* truncated_weight, V* truncated_fraction, V* smallest_ev)
    {
        // Get iteration results from the last sweep
        results_collector& iter = interface_ptr->get_iteration_results();

        *m = 0;
        *truncated_weight = 0;
        *truncated_fraction = 0;
        *smallest_ev = 0;

        if (!iter.empty())
        {
            // iter contains results, one element per microiteration
            const std::vector<boost::any>& m_vec = iter["BondDimension"].get();
            const std::vector<boost::any>& ev_vec = iter["SmallestEV"].get();

            // if we do single-site optimization, we will not have TruncatedWeight or TruncatedFraction, so check whether we have it
            const std::vector<boost::any>& tw_vec = (iter.has("TruncatedWeight")) ? iter["TruncatedWeight"].get() : std::vector<boost::any>();
            const std::vector<boost::any>& tf_vec = (iter.has("TruncatedFraction")) ? iter["TruncatedFraction"].get() : std::vector<boost::any>();

            // return maximum bond dimension
            *m = 0; for (auto&& m_ : m_vec) *m = std::max(*m, boost::any_cast<std::size_t>(m_));
            
            // We return the sum of these values for the last sweep
            // this should be done with transform_reduce
            *truncated_weight = 0; for (auto&& tw_ : tw_vec) *truncated_weight += boost::any_cast<V>(tw_);
            *truncated_fraction = 0; for (auto&& tf_ : tf_vec) *truncated_fraction += boost::any_cast<V>(tf_);
            *smallest_ev = 0; for (auto&& ev_ : ev_vec) *smallest_ev += boost::any_cast<V>(ev_);

            *nsweeps = interface_ptr->get_last_sweep()+1;
        }
        else
            *nsweeps = 0; // If iter is empty, no iterations have been made and thus we return all zeros
    }

    double qcmaquis_interface_get_overlap(const char* filename)
    {
        return interface_ptr->overlap(filename);
    }

    void qcmaquis_interface_stdout(const char* filename)
    {
        stdout_redirect.set_filename(filename);
    }

    void qcmaquis_interface_restore_stdout()
    {
        stdout_redirect.restore();
    }

    int qcmaquis_interface_get_4rdm_elements(int L, const int* slice)
    {
        std::vector<int> slice_ = slice != nullptr ? std::vector<int>(slice, slice+3) : std::vector<int>();
        return measurements_details::get_nrdm_permutations<4>(L, false, slice_);
    }

    int qcmaquis_interface_get_3rdm_elements(int L, bool bra_neq_ket, const int* slice)
    {
        std::vector<int> slice_ = slice != nullptr ? std::vector<int>(slice, slice+1) : std::vector<int>();
        return measurements_details::get_nrdm_permutations<3>(L, bra_neq_ket, slice_);
    }

    void qcmaquis_interface_prepare_hirdm_template(const char* filename, int state, HIRDM_Template tpl, int state_j)
    {
        BaseParameters parms_rdm = parms;
        parms_rdm.erase_measurements();

        // set info for 2U1 group
        #if defined(HAVE_SU2U1PG)
            parms_rdm.set("symmetry", "2u1pg");
        #elif defined(HAVE_SU2U1)
            parms_rdm.set("symmetry", "2u1");
        #endif

        std::string twou1_checkpoint_name;

        int nel = parms_rdm["nelec"];
        int multiplicity = parms_rdm["spin"];

        int Nup, Ndown;

        std::tie(twou1_checkpoint_name, Nup, Ndown) = maquis::interface_detail::twou1_name_Nup_Ndown(pname, state, nel, multiplicity);

        // generate result file name. this should be analogous to SU2U1 file name to be consistent with on-the-fly evaluation
        // which just writes the rdm into the SU2U1 result
        std::string twou1_result_name = std::regex_replace(
                                            maquis::interface_detail::su2u1_name(pname, state), // generate SU2U1 checkpoint file
                                            std::regex("checkpoint"), "results"); // replace 'checkpoint' with 'results'


        // remove the absolute directory for checkpoint and result files
        boost::filesystem::path chkp_name(twou1_checkpoint_name);
        parms_rdm.set("chkpfile", chkp_name.filename().string());

        boost::filesystem::path res_name(twou1_result_name);
        parms_rdm.set("resultfile", res_name.filename().string());

        parms_rdm.set("u1_total_charge1", Nup);
        parms_rdm.set("u1_total_charge2", Ndown);

        if (tpl == TEMPLATE_4RDM)
            parms_rdm.set("MEASURE[4rdm]", "p4:p3:p1:p2@LLL,KKK,III,JJJ");
        else if (tpl == TEMPLATE_TRANSITION_3RDM)
        {
            boost::filesystem::path bra_name(maquis::interface_detail::twou1_name(pname, state_j, nel, multiplicity));
            parms_rdm.set("MEASURE[trans3rdm]", bra_name.filename().string() + ";p1:p2@III,JJJ");
        }
        else
            throw std::runtime_error("Cannot prepare QCMaquis template for this measurement");


        std::ofstream fs(filename);
        fs << parms_rdm;

    }

    void qcmaquis_interface_contract_with_fock_3rdm(double* epsa, int nasht) {
      printf("contract_with_fock epsa = \n");
      for (int i = 0; i < nasht; ++i) {
        printf("%f ", epsa[i]);
      }
      printf("\n");
      // MPS<matrix, TwoU1PG> mps;
      // load(parms["chkpfile"], mps);
      //
      // maquis::integral_map<double> int_map;
      // for (int i = 1; i < nasht + 1; ++i) {
      //   int_map[{i, i, 0, 0}] = epsa[i];
      // }
      // BaseParameters parms_caspt2 = parms;
      // parms.erase("integral_file");
      // parms.erase("integrals");
      // parms.erase("integrals_binary");
      // parms.set("integrals_binary", maquis::serialize(int_map));
      // auto lattice = Lattice(parms_caspt2);
      // auto model = Model<matrix, TwoU1PG>(lattice, parms_caspt2);
      // auto mpo = make_mpo(lattice, model);
      // auto traitClass = MPOTimesMPSTraitClass<tmatrix<double>, TwoU1PG>(
      //     mps, model, lattice, model.total_quantum_numbers(parms),
      //     parms["max_bond_dimension"]);
      // auto outputMPS = traitClass.applyMPO(mpo);
    }
}
