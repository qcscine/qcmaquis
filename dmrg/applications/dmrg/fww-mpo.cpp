/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher
 * Group. See LICENSE.txt for details.
 */

#include "dmrg/utils/DmrgOptions.h"
#include "maquis_dmrg.h"
#include "utils/data_collector.hpp"
#include "utils/io.hpp" // has to be first include because of impi
#include "utils/timings.h"
#include <vector>
#include <cmath>
#include <iostream>
#include <iterator>
#include <memory>
#include <sys/stat.h>
#include <sys/time.h>


#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/sim/interface_sim.h"
#include "dmrg/sim/matrix_types.h"

DmrgParameters build_parms(const maquis::integral_map<double>& int_map) {
  DmrgParameters p;
  p.set("L", 4);
  p.set("symmetry", "2u1pg");
  p.set("u1_total_charge1", 1);
  p.set("u1_total_charge2", 1);
  p.set("nelec", 2);
  p.set("irrep", 0);
  p.set("site_types", "0,0,0,0");
  p.set("integrals_binary", maquis::serialize(int_map));
  p.set("max_bond_dimension", 100);
  p.set("MEASURE[trans3rdm]", "H2.checkpoint_state.0.h5");
  return p;
}

int main(int argc, char **argv) {
  std::cout << "  SCINE QCMaquis \n"
            << "  Quantum Chemical Density Matrix Renormalization group\n"
            << "  available from https://scine.ethz.ch/download/qcmaquis\n"
            << "  based on the ALPS MPS codes from http://alps.comp-phys.org/\n"
            << "  copyright (c) 2015-2018 Laboratory of Physical Chemistry, "
               "ETH Zurich\n"
            << "  copyright (c) 2012-2016 by Sebastian Keller\n"
            << "  copyright (c) 2016-2018 by Alberto Baiardi, Leon Freitag, \n"
            << "  Stefan Knecht, Yingjin Ma \n"
            << "  for details see the publication: \n"
            << "  S. Keller et al., J. Chem. Phys. 143, 244118 (2015)\n"
            << std::endl;

  std::vector<double> fww{-0.71107821123412673, -0.61520367430912049,
                          0.59788654772805638, 0.73289353990937500};
  DmrgOptions opt(argc, argv);
  if (opt.valid) {

    const std::string integrals(
     "    1.63719990472             1     1     1     1\n"
     "  -0.144746632369             1     1     2     1\n"
     "   0.266282775636E-01         2     1     2     1\n"
     "   0.167531821030E-01         2     2     2     1\n"
     "   0.459207160088             1     1     2     2\n"
     "   0.533052674812             2     2     2     2\n"
     "   -0.230673683229E-01        1     1     3     1\n"
     "   0.873924303638E-02         2     1     3     1\n"
     "   0.212121998644E-01         2     2     3     1\n"
     "   0.517422983458E-02         3     1     3     1\n"
     "   0.166602690956E-01         3     2     3     1\n"
     "   0.201358122600E-01         3     3     3     1\n"
     "   0.113600263329             1     1     3     2\n"
     "   0.130338872974E-01         2     1     3     2\n"
     "   0.170776682916             2     2     3     2\n"
     "   0.130799970176             3     2     3     2\n"
     "   0.156301172407             3     3     3     2\n"
     "   0.386363318135             1     1     3     3\n"
     "   0.177101856628E-01         2     1     3     3\n"
     "   0.469384865605             2     2     3     3\n"
     "   0.439202823493             3     3     3     3\n"
     "   -4.96687194130             1     1     0     0\n"
     "   0.128261636279             2     1     0     0\n"
     "   -1.74403215804             2     2     0     0\n"
     "   -0.589618128664E-02        3     1     0     0\n"
     "   -0.377341184513            3     2     0     0\n"
     "   -1.09420984374             3     3     0     0\n"
     "   1.58753163271              0     0     0     0\n");


    DmrgParameters p;
    p.set("integrals",integrals);


    p.set("site_types","0,0,0");
    p.set("L",3);
    p.set("irrep",0);

    p.set("nsweeps",2);
    p.set("max_bond_dimension",100);

    // for SU2U1
    p.set("nelec",4);
    p.set("spin",0);

    // for 2U1
    p.set("symmetry", "2u1pg");
    p.set("u1_total_charge1",2);
    p.set("u1_total_charge2",2);

    // Measure 4-RDM
    // p.set("MEASURE[4rdm]",1);
    // Measure 3-RDM
    // p.set("MEASURE[3rdm]",1);
    
    // checkpoint files
    p.set("chkpfile","H2.checkpoint_state.0.h5");

    std::cout << "Optimizing MPS" << std::endl;
    maquis::DMRGInterface<double> interface(p);
    interface.optimize();

    // auto fourrdm = interface.fourrdm();
    // for (const auto& e : fourrdm.second) {
    //   std::cout << "fourdm = " << e << std::endl;
    // }

    MPS<matrix, TwoU1PG> optMPS;
    load("H2.checkpoint_state.0.h5", optMPS);

    // // double compression_trace;
    // // interface_sim<tmatrix<double>, TwoU1PG> interface(opt.parms);
    // // interface.run("optimize");
    // // auto mps = interface.get_mps();
    // // std::cout << "Input bond dim = " << mps.bond_dimension() << std::endl;
    // // std::cout << mps.description() << std::endl;
    // // std::cout << mps << std::endl;
    // // auto mps1 = compression::l2r_compress(mps, 100, 1e-20, compression_trace, true);
    // // std::cout << "Input COMPRESSED bond dim = " << mps1.bond_dimension() << std::endl;


    // std::cout << "Applying MPO" << std::endl;
    // maquis::integral_map<double> int_map{
    //   {{1, 1, 0, 0}, 1.0}, 
    //   {{2, 2, 0, 0}, 2.0}, 
    //   {{3, 3, 0, 0}, 3.0}, 
    //   {{4, 4, 0, 0}, 4.4}, 
    // };
    DmrgParameters p2(p);

    auto lattice = Lattice(p);
    auto model = Model<matrix, TwoU1PG>(lattice, p);
    auto mpo = make_mpo(lattice, model);
    auto traitClass = MPOTimesMPSTraitClass<tmatrix<double>, TwoU1PG>(
        optMPS, model, lattice, model.total_quantum_numbers(p),
        p["max_bond_dimension"]);
    auto outputMPS = traitClass.applyMPO(mpo);
    std::string MPStimesMPOstr("MPStimesMPO.chkp");
    save(MPStimesMPOstr, outputMPS);
    storage::archive ar(MPStimesMPOstr + "/props.h5", "w");
    ar["/parameters"] << opt.parms;
    ar["/status/sweep"] << 3;
    
    // std::cout << "Measuring" << std::endl;
    p.set("MEASURE[trans3rdm]","MPStimesMPO.chkp");
    p.set("chkpfile", MPStimesMPOstr);
    p.set("resultfile", "results.h5");
    maquis::DMRGInterface<double> interface_measure(p);
    interface_measure.measure();
    
    // std::cout << "Output bond dim = " << outputMPS.bond_dimension() << std::endl;
    // std::cout << outputMPS.description() << std::endl;
    // std::cout << outputMPS << std::endl;
    //
    // auto mps2 = compression::l2r_compress(outputMPS, 100, 1e-20, compression_trace, true);
    // std::cout << "Compressed bond dim = " << mps2.bond_dimension() << std::endl;
    // std::cout << "Compression trace: " << compression_trace << std::endl;
    
    // auto&& meas_trans3rdm = model.measurements().at("transition_threeptdm");

  }
}
