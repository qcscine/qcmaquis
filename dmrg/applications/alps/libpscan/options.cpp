/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "libpscan/options.hpp"
#include <alps/utility/copyright.hpp>
#include <boost/limits.hpp>
#include <boost/throw_exception.hpp>
#include <boost/program_options.hpp>
#include <stdexcept>


Options::Options()
:
time_limit(0.),
use_mpi(false),
valid(false), // shall we really run?
write_xml(false)
{ }

Options::Options(int argc, char** argv)
{
    namespace po = boost::program_options;
    
    programname = std::string(argv[0]);
    valid = true;
    if (argc) {
        std::string filename;
        int dummy;
        po::options_description desc("Allowed options");
        desc.add_options()
        ("help", "produce help message")
        ("license,l", "print license conditions")
        ("mpi", "run in parallel using MPI")
        ("Nmax", po::value<int>(&dummy), "not used, only for backward compatibility with `dmrg`")
        ("time-limit,T", po::value<double>(&time_limit)->default_value(0),"time limit for the simulation")
        ("write-xml","write results to XML files")
        ("input-file", po::value<std::string>(&filename), "input file");
        po::positional_options_description p;
        p.add("input-file", 1);
        
        po::variables_map vm;
        po::store(po::command_line_parser(argc, argv).options(desc).positional(p).run(), vm);
        po::notify(vm);
        
        
        if (vm.count("help")) {
            std::cout << desc << "\n";
            valid=false;
        }
        if (vm.count("license")) {
            alps::print_license(std::cout);
            valid=false;
        }
        if (!valid)
            return;
        
        if (vm.count("mpi")) {
            use_mpi = true;
        }
        
        if (vm.count("write-xml"))
            write_xml = true;
        
        if (!filename.empty())
            jobfilename=boost::filesystem::path(filename);
        else
            boost::throw_exception(std::runtime_error("No job file specified"));
    }
}
