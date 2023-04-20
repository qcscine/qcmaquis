/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "utils/io.hpp"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/version.h"
#include <boost/limits.hpp>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>
#include <stdexcept>


DmrgOptions::DmrgOptions()
:
time_limit(0.),
valid(false)
{ }

DmrgOptions::DmrgOptions(int argc, char** argv)
{
    namespace po = boost::program_options;

    programname = std::string(argv[0]);
    valid = true;
    if (argc) {
        std::string parms_fname, model_fname;

        po::options_description desc("Allowed options");
        desc.add_options()
        ("help,h", "produce help message")
        ("version", "print program version")
        ("time-limit,T", po::value<double>(&time_limit)->default_value(0),"time limit for the simulation")
        ("parms-file", po::value<std::string>(&parms_fname), "input dmrg parameters")
        ("model-file", po::value<std::string>(&model_fname), "input model parameters");
        po::positional_options_description p;
        p.add("parms-file", 1);
        p.add("model-file", 1);

        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);


        if (vm.count("help")) {
            maquis::cout << desc << std::endl;
            maquis::cout << "DMRG Parameters:" << std::endl;
            DmrgParameters().print_description(maquis::cout);
            valid=false;
        }
        if (vm.count("version")) {
            maquis::cout << DMRG_VERSION_STRING << std::endl;
            valid=false;
        }
        if (!valid)
            return;


        /// Load parameters
        std::ifstream param_file(parms_fname.c_str());
        if (!param_file)
            throw std::runtime_error("Could not open parameter file " + parms_fname);
        parms = DmrgParameters(param_file);

        /// Load model parameters from second input (if needed)
        std::string integral_file;
        if (parms.is_set("model_file") && model_fname.empty())
            model_fname = parms["model_file"].str();
        if (!model_fname.empty()) {
            std::ifstream model_ifs(model_fname.c_str());
            if (!model_ifs)
                throw std::runtime_error("Could not open model_parms file.");
            parms << ModelParameters(model_ifs);
        }


        if (!vm["time-limit"].defaulted())
            parms["run_seconds"] = time_limit;

        parallel::params.init();
    }
}
