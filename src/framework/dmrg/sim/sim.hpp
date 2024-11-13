/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#include <boost/algorithm/string.hpp>
#include <boost/algorithm/string/classification.hpp>
#include <boost/algorithm/string/trim.hpp>
#include <cassert>
#include <filesystem>
#include <ostream>
#include <sstream>
#include <string>
#include <vector>

#include "dmrg/block_matrix/symmetry/gsl_coupling.h"
#include "dmrg/block_matrix/symmetry/symmetry_traits.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/measurements.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/sim/sim.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/utils/archive.h"
#include "dmrg/utils/checks.h"
#include "dmrg/utils/parallel/utils.hpp"
#include "dmrg/utils/random.hpp"
#include "dmrg/utils/storage.h"
#include "dmrg/version.h"
#include "utils/io.hpp"

namespace sim_detail {
// Checks if the parameters in the list parm already exists in parms, if not,
// loads it from ar (if it fails, it does nothing)
template <class T>
void load_if_not_exists(
    const T& list, BaseParameters& parms, storage::archive& ar
) {
  BaseParameters tmp_parms;
  ar["/parameters"] >> tmp_parms;
  for (auto&& parm : list) {
    if (!parms.is_set(parm)) {
      if (tmp_parms.is_set(parm)) {
        std::string tmp = tmp_parms[parm];
        parms.set(parm, tmp);
      }
    }
  }
}

template <class SymmGroup>
void print_important_parameters(BaseParameters& parms) {
  bool hasSU2 = symm_traits::HasSU2<SymmGroup>::value;
  bool hasPG = symm_traits::HasPG<SymmGroup>::value;
  bool has2U1 = symm_traits::Has2U1<SymmGroup>::value;

  maquis::cout << "\n === Simulation Parameters ===\n"
               << "     Simulation Type: " << parms["simulation_type"] << '\n'
               << "               Sites: " << parms["L"] << "\n"
               << "            Symmetry: " << parms["symmetry"] << "\n";
  if (hasSU2) {
    maquis::cout << "           Electrons: " << parms["nelec"] << "\n"
                 << "                Spin: " << parms["spin"] << "\n";
  }
  if (has2U1) {
    maquis::cout << "   Spin-up electrons: " << parms["u1_total_charge1"]
                 << "\n"
                 << " Spin-down electrons: " << parms["u1_total_charge2"]
                 << "\n";
  }
  if (hasPG) {
    maquis::cout << "               Irrep: " << parms["irrep"] << std::endl;
  }
  maquis::cout << " =============================\n";
}
}  // namespace sim_detail

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::sim(DmrgParameters& parms_)
    : parms(parms_),
      init_sweep(0),
      init_site(-1),
      restore(false),
      dns((parms["donotsave"] != 0) || !parms.is_set("chkpfile")),
      chkpfile(
          parms.is_set("chkpfile")
              ? boost::trim_right_copy_if(
                    parms["chkpfile"].str(), boost::is_any_of("/ ")
                )
              : ""
      ),
      stop_callback(static_cast<double>(parms["run_seconds"])) {
  maquis::cout << DMRG_VERSION_STRING << std::endl;
  storage::setup(parms);

  bool has2U1 = symm_traits::Has2U1<SymmGroup>::value;
  bool hasPG = symm_traits::HasPG<SymmGroup>::value;
  bool hasSU2 = symm_traits::HasSU2<SymmGroup>::value;

  dmrg_random::engine.seed(parms["seed"]);
  // check possible orbital order in existing MPS before(!) model initialization
  if (!chkpfile.empty()) {
    std::filesystem::path p(chkpfile);
    if (std::filesystem::exists(p) && std::filesystem::exists(p / "props.h5")) {
      maquis::checks::orbital_order_check(parms, chkpfile);
    }
    // Load MPS from checkpoint
    loadMPSAndParams(chkpfile, hasSU2, has2U1, hasPG);
  }
  // TODO(Kalman): This currently does not work for PreBO and Vibrational
  // sim_detail::print_important_parameters<SymmGroup>(parms);

  // Initialise Wigner cache for SU2
  if (hasSU2) {
    initializeWignerCache();
  }

  // Model initialization
  lat = Lattice(parms);
  model = Model<Matrix, SymmGroup>(lat, parms);
  mpo = make_mpo(lat, model);
  all_measurements = model.measurements();
  all_measurements << overlap_measurements<Matrix, SymmGroup>(parms);

  // Final check on the checkpoint MPS after model has been initialised
  // Otherwise, does a fresh MPS initialization
  if (restore) {
    maquis::checks::right_end_check(
        chkpfile, mps, model.total_quantum_numbers(parms)
    );
  } else {
    mps = MPS<Matrix, SymmGroup>(lat.size(), *(model.initializer(lat, parms)));
  }

  all_measurements << autocorrelation_measurements<Matrix, SymmGroup>(
      parms, mps
  );
  assert(mps.length() == lat.size());

  /// Update parameters - after checks have passed
  updateParamsInArchive(chkpfile);

  maquis::cout << "MPS initialization has finished...\n";  // MPS restored now
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::loadMPSAndParams(
    const std::string& chkpfile, const bool hasSU2, const bool has2U1,
    const bool hasPG
) {
  std::filesystem::path p(chkpfile);
  if (std::filesystem::exists(p) && std::filesystem::exists(p / "mps0.h5")) {
    storage::archive ar_in(chkpfile + "/props.h5");
    restore = true;
    if (ar_in.is_scalar("/status/sweep")) {
      ar_in["/status/sweep"] >> init_sweep;

      if (ar_in.is_data("/status/site") && ar_in.is_scalar("/status/site")) {
        ar_in["/status/site"] >> init_site;
      }

      if (init_site == -1) {
        ++init_sweep;
      }

      maquis::cout << "Will start again at site " << init_site << " in sweep "
                   << init_sweep << std::endl;
    }
    // load checkpoint
    maquis::cout << "\n!! WARNING: Checkpoint found. !!\n"
                 << "Continuing calculation from checkpoint " << p.c_str()
                 << "\n\n";
    maquis::checks::symmetry_check(parms, chkpfile);
    load(chkpfile, mps);

    // Try to load some necessary parameters from checkpoint if they're not
    // found in the input file
    if (parms["MODEL"] == "quantum_chemistry") {
      loadParams(ar_in, hasSU2, has2U1, hasPG);
    }
  }
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::loadParams(
    storage::archive& ar_in, const bool hasSU2, const bool has2U1,
    const bool hasPG
) const {
  std::vector<std::string> parms_toload{
      "L", "site_types", "orbital_order", "symmetry"
  };
  if (hasSU2) {
    parms_toload.emplace_back("nelec");
    parms_toload.emplace_back("spin");
  } else if (has2U1) {
    parms_toload.emplace_back("u1_total_charge1");
    parms_toload.emplace_back("u1_total_charge2");
  }
  if (hasPG) {
    parms_toload.emplace_back("irrep");
  }
  // Try loading integrals too, unless integral_file is set
  // TODO: use this also with "integrals"
  if (!parms.is_set("integral_file") && !parms.is_set("integrals")) {
    parms_toload.emplace_back("integrals_binary");
  }
  //
  sim_detail::load_if_not_exists(parms_toload, parms, ar_in);
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::initializeWignerCache() const {
  if (!(parms.is_set("NoWignerCache") && parms["NoWignerCache"])) {
    WignerWrapper::UseCache = true;
    int nelec = parms["nelec"];
    int L = parms["L"];
    int spin = parms["spin"];
    int max_spin = (nelec > L) ? nelec - (nelec - L) * 2 : nelec;
    // maximum parameter for 9j symbols
    int max_9j = (max_spin + spin) / 2;
    // For spin = 0 or 1 we may require 9j symbols with c,g,f=2 due to
    // permutational symmetry of the rows/columns since we do not implement this
    // permutational symmetry for performance reasons, we need to fill the cache
    // for elements up to 2
    if (max_9j < 2) {
      max_9j = 2;
    }
    WignerWrapper::fill_cache(max_9j);
  }
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::updateParamsInArchive(const std::string& chkpfile
) const {
  if (!rfile().empty()) {
    storage::archive ar(rfile(), "w");
    ar["/parameters"] << parms;
    ar["/version"] << DMRG_VERSION_STRING;
  }
  if (!dns && !chkpfile.empty()) {
    if (!std::filesystem::exists(chkpfile)) {
      std::filesystem::create_directory(chkpfile);
    }
    storage::archive ar(chkpfile + "/props.h5", "w");

    ar["/parameters"] << parms;
    ar["/version"] << DMRG_VERSION_STRING;
  }
}

template <class Matrix, class SymmGroup>
typename sim<Matrix, SymmGroup>::measurements_type
sim<Matrix, SymmGroup>::iteration_measurements(int sweep) {
  measurements_type mymeas(all_measurements);
  mymeas << overlap_measurements<Matrix, SymmGroup>(parms, sweep);

  measurements_type sweep_measurements;
  if (!parms["ALWAYS_MEASURE"].empty()) {
    sweep_measurements = meas_sublist(mymeas, parms["ALWAYS_MEASURE"]);
  }

  return sweep_measurements;
}

template <class Matrix, class SymmGroup>
sim<Matrix, SymmGroup>::~sim() = default;

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::checkpoint_simulation(
    const MPS<Matrix, SymmGroup>& state, const status_type& status,
    std::string filename
) {
  std::string chkpfilename;
  if (filename.empty()) {
    chkpfilename = chkpfolder();
  } else {
    chkpfilename = chkpfolder() + "_" + filename;
  }
  if (!dns && !chkpfilename.empty()) {
    /// save state to chkp dir
    save(chkpfilename, state);

    /// save status
    if (!parallel::master()) {
      return;
    }
    storage::archive ar(chkpfilename + "/props.h5", "w");
    ar["/status"] << status;
  }
}

template <class Matrix, class SymmGroup>
std::string sim<Matrix, SymmGroup>::results_archive_path(
    const status_type& status
) const {
  std::ostringstream oss;
  oss.str("");
#if defined(__xlC__) || defined(__FCC_VERSION)
  typename status_type::const_iterator match = status.find("sweep");
  oss << "/spectrum/iteration/" << match->second;
#else
  oss << "/spectrum/iteration/" << status.at("sweep");
#endif
  return oss.str();
}

template <class Matrix, class SymmGroup>
void sim<Matrix, SymmGroup>::measure(
    std::string archive_path, measurements_type& meas
) {
#ifdef MAQUIS_OPENMP
  if (parms["parallelize_measurements"]) {
#pragma omp parallel for schedule(dynamic)
    for (typename measurements_type::iterator it = meas.begin();
         it < meas.end(); it++) {
      MPS<Matrix, SymmGroup> mpsCopy =
          mps;  // this is required as the measurements might change the pairing
                // of the mps
      // note that omp firstprivate cannot be used since the mps does apparently
      // not fulfill the necessary requirements
      measure_and_save<Matrix, SymmGroup> ms(rfile(), archive_path, mpsCopy);
      ms(*it);
    }
  } else
    std::for_each(
        meas.begin(), meas.end(),
        measure_and_save<Matrix, SymmGroup>(rfile(), archive_path, mps)
    );
#else
  std::for_each(
      meas.begin(), meas.end(),
      measure_and_save<Matrix, SymmGroup>(rfile(), archive_path, mps)
  );
#endif

  // TODO: move into special measurement
  std::vector<int>* measure_es_where = nullptr;
  entanglement_spectrum_type* spectra = nullptr;
  if (parms.defined("entanglement_spectra")) {
    spectra = new entanglement_spectrum_type();
    measure_es_where = new std::vector<int>();
    *measure_es_where =
        parms.template get<std::vector<int> >("entanglement_spectra");
  }
  std::vector<double> entropies;
  std::vector<double> renyi2;
  if (parms["MEASURE[Entropy]"]) {
    maquis::cout << "Calculating vN entropy." << std::endl;
    entropies = calculate_bond_entropies(mps);
  }
  if (parms["MEASURE[Renyi2]"]) {
    maquis::cout << "Calculating n=2 Renyi entropy." << std::endl;
    renyi2 = calculate_bond_renyi_entropies(mps, 2, measure_es_where, spectra);
  }

  if (!rfile().empty()) {
    storage::archive ar(rfile(), "w");
    if (!entropies.empty()) {
      ar[archive_path + "Entropy/mean/value"] << entropies;
    }
    if (!renyi2.empty()) {
      ar[archive_path + "Renyi2/mean/value"] << renyi2;
    }
    if (spectra != nullptr) {
      ar[archive_path + "Entanglement Spectra/mean/value"] << *spectra;
    }
  }
}
