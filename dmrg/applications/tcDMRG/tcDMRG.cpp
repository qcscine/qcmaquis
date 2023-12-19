/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2023 Institute for Theoretical Physics, ETH Zurich
 *               2023- by Alberto Baiardi <lefreita@ethz.ch>
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
#include <iostream>
#include <sys/stat.h>
#include <sys/time.h>
#include "utils/data_collector.hpp"
#include "utils/timings.h"
#include "maquis_dmrg.h"
#include "dmrg/utils/DmrgOptions.h"


int main(int argc, char ** argv)
{
  std::cout << "  SCINE QCMaquis " << std::endl
            << "  Quantum Chemical Density Matrix Renormalization group" << std::endl
            << "  available from https://scine.ethz.ch/download/qcmaquis" << std::endl
            << "  based on the ALPS MPS codes from http://alps.comp-phys.org" << std::endl
            << "  copyright (c) 2015-2018 Laboratory of Physical Chemistry, ETH Zurich" << std::endl
            << "  copyright (c) 2012-2016 by Sebastian Keller" << std::endl
            << "  copyright (c) 2016-2021 by Alberto Baiardi, Leon Freitag," << std::endl
            << "  Stefan Knecht, Yingjin Ma" << std::endl
            << "  For details on tcDMRG see the following publications:" << std::endl
            << "   - A. Baiardi, M. Reiher, J. Chem. Phys. 153, 164115 (2020)" << std::endl
            << "   - A. Baiardi, M. Lesiuk, M. Reiher, J. Chem. Theory Comput. 18, 4203 (2022)" << std::endl;
  DmrgOptions opt(argc, argv);
  if (opt.valid) {
    maquis::cout.precision(10);
    DCOLLECTOR_SET_SIZE(gemm_collector, opt.parms["max_bond_dimension"]+1)
    DCOLLECTOR_SET_SIZE(svd_collector, opt.parms["max_bond_dimension"]+1)
    timeval now, then, snow, sthen;
    gettimeofday(&now, NULL);
    maquis::DMRGInterface<double> interface(opt.parms);
    interface.runTranscorrelated();
    gettimeofday(&then, NULL);
    double elapsed = then.tv_sec-now.tv_sec + 1e-6 * (then.tv_usec-now.tv_usec);
    DCOLLECTOR_SAVE_TO_FILE(gemm_collector, "collectors.h5", "/results")
    DCOLLECTOR_SAVE_TO_FILE(svd_collector, "collectors.h5", "/results")
    maquis::cout << "Task took " << elapsed << " seconds." << std::endl;
  }
}