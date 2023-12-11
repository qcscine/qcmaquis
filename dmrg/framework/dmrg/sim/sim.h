/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef APP_SIM_H
#define APP_SIM_H

#include <cmath>
#include <iterator>
#include <iostream>
#include "utils/data_collector.hpp"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/twositetensor.h"
#include "dmrg/mp_tensors/contractions.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/utils/random.hpp"
#include "dmrg/utils/time_stopper.h"
#include "utils/timings.h"
#include "dmrg/utils/checks.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "abstract_sim.h"

template <class Matrix, class SymmGroup>
class sim : public abstract_sim {
private:
  void initializeWignerCache() const;
  void loadParams(
      storage::archive& ar_in,
      const bool hasSU2, const bool has2U1, const bool hasPG) const;
  void loadMPSAndParams(
      const std::string& chkpfile,
      const bool hasSU2, const bool has2U1, const bool hasPG);
  void updateParamsInArchive(const std::string& chkpfile) const;
public:
    explicit sim(DmrgParameters &);
     ~sim() override;

protected:
    using measurements_type = typename Model<Matrix, SymmGroup>::measurements_type;
    using status_type = std::map<std::string, int>;

    virtual std::string results_archive_path(status_type const&) const;

    measurements_type iteration_measurements(int sweep);
    virtual void measure(std::string archive_path, measurements_type & meas);
    // TODO: can be made const, now only problem are parameters
    virtual void checkpoint_simulation(MPS<Matrix, SymmGroup> const& state, status_type const&, std::string filename = "");

    DmrgParameters& parms;
    int init_sweep, init_site;
    bool restore;
    bool dns;
    std::string chkpfile;
    std::string chkpfolder() const {
      if (parms.is_set("chkpfile")) {
        return parms["chkpfile"].str();
      } else {
        return std::string();
}
    }
    std::string rfile() const
    {
      if (parms.is_set("resultfile")) {
        return parms["resultfile"].str();
      } else {
        return std::string();
      }
    };
    time_stopper stop_callback;
    Lattice lat;
    Model<Matrix, SymmGroup> model;
    MPS<Matrix, SymmGroup> mps;
    MPO<Matrix, SymmGroup> mpo, mpoc;
    measurements_type all_measurements, sweep_measurements;
};


#include "sim.hpp"
#endif
