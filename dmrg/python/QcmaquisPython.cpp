#include "dmrg/tools/determinant.h"
#include "dmrg/tools/electronic_srcas.h"
#include "dmrg/tools/vib_srcas.h"
#include "integral_helper.h"
#include "maquis_dmrg.h"
#include <pybind11/complex.h>
#include <pybind11/detail/common.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <complex>

namespace py = pybind11;

// TODO: bind results_collector meas_with_results_type
PYBIND11_MODULE(_dmrg, module) {
  using ComplexElectronicTranscorrelatedIntegralMap =
      chem::integral_map<std::complex<double>, chem::Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>;
  using RealElectronicTranscorrelatedIntegralMap =
      chem::integral_map<double, chem::Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>;

  module.def("asymmetry",
             [](const std::string& str) { // Accepts str or bytes from Python
               std::cout << str << "\n";
               std::cout << (str.c_str()) << "\n";
               return str.c_str(); // Looks harmless, but implicitly converts to str
             });

  /* Transcorrelated Complex Integral Map */
  py::class_<ComplexElectronicTranscorrelatedIntegralMap>(module, "ComplexTCIntegralMap")
      .def(py::init<>())
      .def("set", &ComplexElectronicTranscorrelatedIntegralMap::set)
      .def("at", [](ComplexElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah,
                    double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const ComplexElectronicTranscorrelatedIntegralMap& blub,
                    const std::array<int, 6>& blah) { return &blub.at(blah); })
      .def("__getitem__",
           [](ComplexElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah) { return &blub[blah]; })
      .def("__setitem__",
           [](ComplexElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah) { return &blub[blah]; });

  /* Transcorrelated Real Integral Map */
  py::class_<RealElectronicTranscorrelatedIntegralMap>(module, "TCIntegralMap")
      .def(py::init<>())
      .def("set", &RealElectronicTranscorrelatedIntegralMap::set)
      .def("at", [](RealElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah,
                    double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const RealElectronicTranscorrelatedIntegralMap& blub,
                    const std::array<int, 6>& blah) { return &blub.at(blah); })
      .def("__getitem__",
           [](RealElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah) { return &blub[blah]; })
      .def("__setitem__",
           [](RealElectronicTranscorrelatedIntegralMap& blub, const std::array<int, 6>& blah) { return &blub[blah]; });

  /* Real Integral Map */
  py::class_<chem::integral_map<double>>(module, "IntegralMap")
      .def(py::init<>())
      .def("set", &chem::integral_map<double>::set)
      .def("at",
           [](chem::integral_map<double>& blub, const std::array<int, 4>& blah, double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub.at(blah); })
      .def("__getitem__", [](chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub[blah]; })
      .def("__setitem__", [](chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub[blah]; });

  /* Complex Integral Map */
  py::class_<chem::integral_map<std::complex<double>>>(module, "ComplexIntegralMap")
      .def(py::init<>())
      .def("set", &chem::integral_map<std::complex<double>>::set)
      .def("at", [](chem::integral_map<std::complex<double>>& blub, const std::array<int, 4>& blah,
                    std::complex<double>& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const chem::integral_map<std::complex<double>>& blub,
                    const std::array<int, 4>& blah) { return &blub.at(blah); })
      .def("__getitem__",
           [](chem::integral_map<std::complex<double>>& blub, const std::array<int, 4>& blah) { return &blub[blah]; })
      .def("__setitem__",
           [](chem::integral_map<std::complex<double>>& blub, const std::array<int, 4>& blah) { return &blub[blah]; });

  /* BaseParameters for DmrgParameters */
  py::class_<BaseParameters>(module, "BaseParameters")
      .def(py::init<>())
      .def("set", &BaseParameters::set<int>)
      .def("set", &BaseParameters::set<std::string>)
      .def("set", &BaseParameters::set<double>)
      .def("erase", &BaseParameters::erase);
  /* DmrgParameters */
  py::class_<DmrgParameters, BaseParameters>(module, "DmrgParameters").def(py::init<>());

  /* Base Real Electronic SRCAS class */
  py::class_<maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>>(module, "BaseRealElectronicSRCAS")
      .def("run", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::run)
      .def("currentQueen", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::currentQueen)
      .def("sampledTable", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::sampledTable)
      .def("completeness", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::completeness)
      .def("printSettings", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::printSettings)
      .def("printResults", &maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>::printResults);
  /* Real Electronic SRCAS */
  py::class_<maquis::srcas::ElectronicSRCAS<double>, maquis::srcas::BaseSRCAS<double, maquis::srcas::Determinant>>(
      module, "RealElectronicSRCAS")
      .def(py::init<DmrgParameters&, std::shared_ptr<maquis::DMRGInterface<double>>>());

  /* Base Complex Electronic SRCAS class */
  py::class_<maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>>(module,
                                                                                         "BaseComplexElectronicSRCAS")
      .def("run", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::run)
      .def("currentQueen", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::currentQueen)
      .def("sampledTable", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::sampledTable)
      .def("completeness", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::completeness)
      .def("printSettings", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::printSettings)
      .def("printResults", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>::printResults);
  /* Complex Electronic SRCAS */
  py::class_<maquis::srcas::ElectronicSRCAS<std::complex<double>>, maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::Determinant>>(
      module, "ComplexElectronicSRCAS")
      .def(py::init<DmrgParameters&, std::shared_ptr<maquis::DMRGInterface<std::complex<double>>>>());

  /* Base Real Electronic SRCAS class */
  py::class_<maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>>(module, "BaseRealVibSRCAS")
      .def("run", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::run)
      .def("currentQueen", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::currentQueen)
      .def("sampledTable", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::sampledTable)
      .def("completeness", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::completeness)
      .def("printSettings", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::printSettings)
      .def("printResults", &maquis::srcas::BaseSRCAS<double, maquis::srcas::VibONV>::printResults);
  /* Real Vibrational SRCAS */
  py::class_<maquis::srcas::VibSRCAS<double>>(module, "RealVibrationalSRCAS")
      .def(py::init<DmrgParameters&, std::shared_ptr<maquis::DMRGInterface<double>>>());

  /* Base Complex Electronic SRCAS class */
  py::class_<maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>>(module, "BaseComplexVibSRCAS")
      .def("run", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::run)
      .def("currentQueen", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::currentQueen)
      .def("sampledTable", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::sampledTable)
      .def("completeness", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::completeness)
      .def("printSettings", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::printSettings)
      .def("printResults", &maquis::srcas::BaseSRCAS<std::complex<double>, maquis::srcas::VibONV>::printResults);
  /* Complex Vibrational SRCAS */
  py::class_<maquis::srcas::VibSRCAS<std::complex<double>>>(module, "ComplexVibrationalSRCAS")
      .def(py::init<DmrgParameters&, std::shared_ptr<maquis::DMRGInterface<std::complex<double>>>>());

  /* Real Electronic DMRG Interface */
  py::class_<maquis::DMRGInterface<double>, std::shared_ptr<maquis::DMRGInterface<double>>>(module, "DmrgReal")
      .def(py::init<DmrgParameters&>())
      .def("optimize", &maquis::DMRGInterface<double>::optimize)
      .def("evolve", &maquis::DMRGInterface<double>::evolve)
      .def("runInversePowerIteration", &maquis::DMRGInterface<double>::runInversePowerIteration)
      .def("runFEAST", &maquis::DMRGInterface<double>::runFEAST)
      .def("energy", &maquis::DMRGInterface<double>::energy)
      .def("energyFEAST", &maquis::DMRGInterface<double>::energyFEAST)
      .def("getCICoefficient", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::getCICoefficient)
      .def("get_iteration_results", &maquis::DMRGInterface<double>::get_iteration_results)
      .def("get_last_sweep", &maquis::DMRGInterface<double>::get_last_sweep)
      .def("run_measure", &maquis::DMRGInterface<double>::run_measure)
      .def("measure", &maquis::DMRGInterface<double>::measure)
      .def("measurements", &maquis::DMRGInterface<double>::measurements)
      .def("update_tc_integrals", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::update_tc_integrals)
      .def("update_integrals", py::overload_cast<std::string>(&maquis::DMRGInterface<double>::update_integrals))
      .def("update_integrals",
           py::overload_cast<const chem::integral_map<double>&>(&maquis::DMRGInterface<double>::update_integrals))
      .def("fiedler_order", &maquis::DMRGInterface<double>::fiedler_order)
      .def("onerdm", &maquis::DMRGInterface<double>::onerdm)
      .def("onespdm", &maquis::DMRGInterface<double>::onespdm)
      .def("twordm", &maquis::DMRGInterface<double>::twordm)
      .def("threerdm", &maquis::DMRGInterface<double>::threerdm)
      .def("fourrdm", &maquis::DMRGInterface<double>::fourrdm)
      .def("getMeasurement", &maquis::DMRGInterface<double>::getMeasurement)
      .def("measure_and_save_3drm", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::measure_and_save_3rdm)
      .def("measure_and_save_4drm", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::measure_and_save_4rdm)
      .def("mutinf", &maquis::DMRGInterface<double>::mutinf)
      .def("measure_and_save_trans3rdm", &maquis::DMRGInterface<double>::measure_and_save_trans3rdm)
      .def("overlap", &maquis::DMRGInterface<double>::overlap)
      .def("dump_parameters", &maquis::DMRGInterface<double>::dump_parameters);

  /* Complex Electronic DMRG Interface */
  py::class_<maquis::DMRGInterface<std::complex<double>>, std::shared_ptr<maquis::DMRGInterface<std::complex<double>>>>(
      module, "DmrgComplex")
      .def(py::init<DmrgParameters&>())
      .def("optimize", &maquis::DMRGInterface<std::complex<double>>::optimize)
      .def("evolve", &maquis::DMRGInterface<std::complex<double>>::evolve)
      .def("runInversePowerIteration", &maquis::DMRGInterface<std::complex<double>>::runInversePowerIteration)
      .def("runFEAST", &maquis::DMRGInterface<std::complex<double>>::runFEAST)
      .def("energy", &maquis::DMRGInterface<std::complex<double>>::energy)
      .def("energyFEAST", &maquis::DMRGInterface<std::complex<double>>::energyFEAST)
      .def("getCICoefficient", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::getCICoefficient)
      .def("get_iteration_results", &maquis::DMRGInterface<std::complex<double>>::get_iteration_results)
      .def("get_last_sweep", &maquis::DMRGInterface<std::complex<double>>::get_last_sweep)
      .def("run_measure", &maquis::DMRGInterface<std::complex<double>>::run_measure)
      .def("measure", &maquis::DMRGInterface<std::complex<double>>::measure)
      .def("measurements", &maquis::DMRGInterface<std::complex<double>>::measurements)
      .def("update_tc_integrals", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::update_tc_integrals)
      .def("update_integrals", py::overload_cast<std::string>(&maquis::DMRGInterface<std::complex<double>>::update_integrals))
      .def("update_integrals", py::overload_cast<const chem::integral_map<std::complex<double>>&>(
                                   &maquis::DMRGInterface<std::complex<double>>::update_integrals))
      .def("fiedler_order", &maquis::DMRGInterface<std::complex<double>>::fiedler_order)
      .def("onerdm", &maquis::DMRGInterface<std::complex<double>>::onerdm)
      .def("onespdm", &maquis::DMRGInterface<std::complex<double>>::onespdm)
      .def("twordm", &maquis::DMRGInterface<std::complex<double>>::twordm)
      .def("threerdm", &maquis::DMRGInterface<std::complex<double>>::threerdm)
      .def("fourrdm", &maquis::DMRGInterface<std::complex<double>>::fourrdm)
      .def("getMeasurement", &maquis::DMRGInterface<std::complex<double>>::getMeasurement)
      .def("measure_and_save_3drm",
           &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::measure_and_save_3rdm)
      .def("measure_and_save_4drm",
           &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::measure_and_save_4rdm)
      .def("mutinf", &maquis::DMRGInterface<std::complex<double>>::mutinf)
      .def("measure_and_save_trans3rdm", &maquis::DMRGInterface<std::complex<double>>::measure_and_save_trans3rdm)
      .def("overlap", &maquis::DMRGInterface<std::complex<double>>::overlap)
      .def("dump_parameters", &maquis::DMRGInterface<std::complex<double>>::dump_parameters);
}
