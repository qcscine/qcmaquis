/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "utils/io.hpp" // has to be first include because of impi
#include "dmrg/utils/DmrgOptions.h"
#include "simulation.h"
#include "dmrg/sim/symmetry_factory.h"

#include <iostream>
#include <sys/time.h>
#include <sys/stat.h>

#include <mpi_interface.h>
#include <utils/qcmaquis_header.h>

namespace maquis
{
  Scine::Mpi::MpiInterface* mpi__;
  std::unique_ptr<Scine::Mpi::MpiInterface> mpi;
}


int main(int argc, char ** argv)
{

    // setup MPI interface. It does nothing for serial runs
    if (!maquis::mpi__) {
        maquis::mpi   = std::unique_ptr<Scine::Mpi::MpiInterface>(new Scine::Mpi::MpiInterface(nullptr, 0));
        maquis::mpi__ = maquis::mpi.get();
    }


    if(maquis::mpi__->getGlobalRank() == 0)
        maquis::qcmaquis_header(maquis::mpi__->isMPIAvailable(),maquis::mpi__->getGlobalCommunicatorSize());

    DmrgOptions opt(argc, argv);
    if (opt.valid) {

        maquis::cout.precision(10);

        timeval now, then, snow, sthen;
        gettimeofday(&now, NULL);


        try {
            simulation_traits::shared_ptr sim = dmrg::symmetry_factory<simulation_traits>(opt.parms);
            sim->run(opt.parms);
        } catch (std::exception & e) {
            maquis::cerr << "Exception thrown!" << std::endl;
            maquis::cerr << e.what() << std::endl;
            exit(1);
        }

        gettimeofday(&then, NULL);
        double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);

        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
    }

    // terminate MPI (does nothing if serial run)
    maquis::mpi.reset(nullptr);
    maquis::mpi__ = nullptr;

}
