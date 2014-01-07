/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/block_matrix/symmetry/nu1.h"
typedef NU1 grp;

#include "dmrg/block_matrix/detail/alps.hpp"
typedef alps::numeric::matrix<double>                matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

#include "tevol_sim.hpp"
#include "dmrg/evolve/tevol_nn_sim.h"
#include "dmrg/evolve/tevol_mpo_sim.h"

#include <boost/filesystem/fstream.hpp>
#include <string>


#include "libpscan/run_sim.hpp"

void run_sim(const boost::filesystem::path& infile, const boost::filesystem::path& outfile,
             double time_limit)
{
    maquis::cout.precision(10);
    
    /// Load parameters
    DmrgParameters parms;
    {
        boost::filesystem::ifstream param_file(infile);
        if (!param_file)
            throw std::runtime_error(std::string("Could not open parameter file ") + infile.string() +".");
        alps::Parameters p; p.extract_from_xml(param_file);
        parms = DmrgParameters(p);
    }
    
    /// Match parameters of ALPS DMRG
    parms.set("nsweeps", int(parms["SWEEPS"]));
    parms.set("max_bond_dimension", int(parms["MAXSTATES"]));
    parms.set("chkpfile",   outfile.stem().string() + ".chkp");
    parms.set("resultfile", outfile.stem().string() + ".h5");
    parms.set("run_seconds", time_limit);
    
    /// Check which matrix to use and which time evolution
    boost::scoped_ptr<abstract_sim> sim;
    if (parms["COMPLEX"]) {
        if (parms["te_type"] == "nn")
            sim.reset(new tevol_sim<cmatrix, grp, nearest_neighbors_evolver<cmatrix, grp> >(parms));
        else if (parms["te_type"] == "mpo")
            sim.reset(new tevol_sim<cmatrix, grp, mpo_evolver<cmatrix, grp> >(parms));
    } else {
        if (parms["te_type"] == "nn")
            sim.reset(new tevol_sim<matrix, grp, nearest_neighbors_evolver<matrix, grp> >(parms));
        else if (parms["te_type"] == "mpo")
            sim.reset(new tevol_sim<matrix, grp, mpo_evolver<matrix, grp> >(parms));
    }
    
    /// Run
    sim->run();
}
