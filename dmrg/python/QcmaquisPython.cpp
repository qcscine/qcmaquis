#include "integral_helper.h"
#include "maquis_dmrg.h"
#include <pybind11/complex.h>
#include <pybind11/detail/common.h>
#include <pybind11/iostream.h>
#include <pybind11/pybind11.h>
#include <pybind11/stl.h>
#include <complex>
// #include <pybind11/stl_bind.h>

namespace py = pybind11;

// TODO: bind results_collector meas_with_results_type
PYBIND11_MODULE(_dmrg, module) {
  // py::bind_map<std::unordered_map<std::pair<std::string, std::string>, double, boost::hash<std::pair<std::string,
  // std::string>>>>(
  //     m, "CIMap");
  module.def("asymmetry",
             [](std::string str) { // Accepts str or bytes from Python
               std::cout << str << std::endl;
               std::cout << (str.c_str()) << std::endl;
               return str.c_str(); // Looks harmless, but implicitly converts to str
             });
  py::class_<chem::integral_map<std::complex<double>, chem::Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>>(module, "ComplexTCIntegralMap")
      .def(py::init<>())
      .def("set", &chem::integral_map<std::complex<double>, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>::set)
      .def("at", [](chem::integral_map<std::complex<double>, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah, double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const chem::integral_map<std::complex<double>, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub.at(blah); })
      .def("__getitem__", [](chem::integral_map<std::complex<double>, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub[blah]; })
      .def("__setitem__", [](chem::integral_map<std::complex<double>, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub[blah]; });


  py::class_<chem::integral_map<double, chem::Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>>(module, "TCIntegralMap")
      .def(py::init<>())
      .def("set", &chem::integral_map<double, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>::set)
      .def("at", [](chem::integral_map<double, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah, double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const chem::integral_map<double, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub.at(blah); })
      .def("__getitem__", [](chem::integral_map<double, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub[blah]; })
      .def("__setitem__", [](chem::integral_map<double, Hamiltonian::Electronic, chem::HamiltonianTransformation::Transcorrelated>& blub, const std::array<int, 6>& blah) { return &blub[blah]; });

  // integral_map
  py::class_<chem::integral_map<double>>(module, "IntegralMap")
      .def(py::init<>())
      .def("set", &chem::integral_map<double>::set)
      .def("at", [](chem::integral_map<double>& blub, const std::array<int, 4>& blah, double& hihi) { blub.at(blah) = hihi; })
      .def("at", [](const chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub.at(blah); })
      .def("__getitem__", [](chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub[blah]; })
      .def("__setitem__", [](chem::integral_map<double>& blub, const std::array<int, 4>& blah) { return &blub[blah]; });

  // BaseParameters for DmrgParameters
  py::class_<BaseParameters>(module, "BaseParameters")
      .def(py::init<>())
      .def("set", &BaseParameters::set<int>)
      .def("set", &BaseParameters::set<std::string>)
      .def("set", &BaseParameters::set<double>)
      .def("erase", &BaseParameters::erase);

  // DmrgParameters
  py::class_<DmrgParameters, BaseParameters>(module, "DmrgParameters").def(py::init<>());

  // Dmrg real
  py::class_<maquis::DMRGInterface<double>>(module, "DmrgReal")
      .def(py::init<DmrgParameters&>())
      .def("optimize", &maquis::DMRGInterface<double>::optimize)
      .def("evolve", &maquis::DMRGInterface<double>::evolve)
           // py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
      .def("energy", &maquis::DMRGInterface<double>::energy)
      .def("get_iteration_results", &maquis::DMRGInterface<double>::get_iteration_results)
      .def("get_last_sweep", &maquis::DMRGInterface<double>::get_last_sweep)
      .def("run_measure", &maquis::DMRGInterface<double>::run_measure)
      .def("measure", &maquis::DMRGInterface<double>::measure)
      .def("measurements", &maquis::DMRGInterface<double>::measurements)
      .def("update_integrals", py::overload_cast<std::string>(&maquis::DMRGInterface<double>::update_integrals))
      .def("update_integrals", py::overload_cast<const chem::integral_map<double>&>(&maquis::DMRGInterface<double>::update_integrals))
      .def("update_tc_integrals", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::update_tc_integrals)
      // .def("fiedler_order", &maquis::DMRGInterface<double>::fiedler_order)
      .def("onerdm", &maquis::DMRGInterface<double>::onerdm)
      .def("onespdm", &maquis::DMRGInterface<double>::onespdm)
      .def("twordm", &maquis::DMRGInterface<double>::twordm)
      .def("threerdm", &maquis::DMRGInterface<double>::threerdm)
      .def("fourrdm", &maquis::DMRGInterface<double>::fourrdm)
      .def("measure_and_save_3drm", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::measure_and_save_3rdm)
      .def("measure_and_save_4drm", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::measure_and_save_4rdm)
      .def("mutinf", &maquis::DMRGInterface<double>::mutinf)
      .def("measure_and_save_trans3rdm", &maquis::DMRGInterface<double>::measure_and_save_trans3rdm)
      .def("overlap", &maquis::DMRGInterface<double>::overlap)
      // .def("getCICoefficients", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::getCICoefficients)
      .def("getCICoefficient", &maquis::DMRGInterface<double, chem::Hamiltonian::Electronic>::getCICoefficient)
      .def("dump_parameters", &maquis::DMRGInterface<double>::dump_parameters);

  // Dmrg complex
  py::class_<maquis::DMRGInterface<std::complex<double>>>(module, "DmrgComplex")
      .def(py::init<DmrgParameters&>())
      .def("optimize", &maquis::DMRGInterface<std::complex<double>>::optimize)
      // .def("runFEAST", &maquis::DMRGInterface<std::complex<double>>::runFEAST)
      .def("evolve", &maquis::DMRGInterface<std::complex<double>>::evolve)
           // py::call_guard<py::scoped_ostream_redirect, py::scoped_estream_redirect>())
      .def("energy", &maquis::DMRGInterface<std::complex<double>>::energy)
      // .def("energyFEAST", &maquis::DMRGInterface<std::complex<double>>::energyFEAST)
      .def("get_iteration_results", &maquis::DMRGInterface<std::complex<double>>::get_iteration_results)
      .def("get_last_sweep", &maquis::DMRGInterface<std::complex<double>>::get_last_sweep)
      .def("run_measure", &maquis::DMRGInterface<std::complex<double>>::run_measure)
      .def("measure", &maquis::DMRGInterface<std::complex<double>>::measure)
      .def("measurements", &maquis::DMRGInterface<std::complex<double>>::measurements)
      .def("update_integrals", py::overload_cast<std::string>(&maquis::DMRGInterface<std::complex<double>>::update_integrals))
      .def("update_integrals", py::overload_cast<const chem::integral_map<std::complex<double>>&>(&maquis::DMRGInterface<std::complex<double>>::update_integrals))
      .def("update_tc_integrals", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::update_tc_integrals)
      // .def("fiedler_order", &maquis::DMRGInterface<std::complex<double>>::fiedler_order)
      .def("onerdm", &maquis::DMRGInterface<std::complex<double>>::onerdm)
      .def("onespdm", &maquis::DMRGInterface<std::complex<double>>::onespdm)
      .def("twordm", &maquis::DMRGInterface<std::complex<double>>::twordm)
      .def("threerdm", &maquis::DMRGInterface<std::complex<double>>::threerdm)
      .def("fourrdm", &maquis::DMRGInterface<std::complex<double>>::fourrdm)
      // .def("measure_and_save_3drm", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::measure_and_save_3rdm)
      // .def("measure_and_save_4drm", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::measure_and_save_4rdm)
      .def("mutinf", &maquis::DMRGInterface<std::complex<double>>::mutinf)
      .def("measure_and_save_trans3rdm", &maquis::DMRGInterface<std::complex<double>>::measure_and_save_trans3rdm)
      .def("overlap", &maquis::DMRGInterface<std::complex<double>>::overlap)
      // .def("getCICoefficients", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::getCICoefficients)
      .def("getCICoefficient", &maquis::DMRGInterface<std::complex<double>, chem::Hamiltonian::Electronic>::getCICoefficient)
      .def("dump_parameters", &maquis::DMRGInterface<std::complex<double>>::dump_parameters);

  // py::class_<maquis::FiedlerGenerator<double>>(module, "FiedlerGenerator")
  //     .def(py::init<DmrgParameters&, std::vector<std::vector<int>>, int>())
  //     .def("get_fiedler_order", &maquis::FiedlerGenerator<double>::getFiedlerOrder);
}
