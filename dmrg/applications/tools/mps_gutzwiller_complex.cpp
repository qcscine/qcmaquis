/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#include <iostream>
#include <string>
#include "GutzwillerCorrelatorClass.h"
#include "dmrg/sim/matrix_types.h"

// Decides the symmetry type
#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#endif

int main(int argc, char** argv) {
  // Info regarding the code usage
  if (argc != 3) {
    maquis::cout << "Usage: ./mps_apply_Gutzwiller <mps.h5> <JValue>"
                 << std::endl;
    maquis::cout
        << "Applies the Gutzwiller correlator onto the MPS wave function "
        << std::endl;
    maquis::cout << "See J. Chem. Phys. 153, 164115 (2020)" << std::endl;
    exit(1);
  }
  // Actual calculation
  maquis::cout.precision(10);
  std::ifstream param_file(argv[1]);
  if (!param_file) {
    maquis::cerr << "Could not open the mps." << std::endl;
  } else {
    auto jValue = std::stod(argv[2]);
    GutzwillerCalculator<cmatrix, grp> calculator(argv[1], jValue);
    calculator.applyCorrelator();
    auto mpsApplied = calculator.getMPS();
    // Save the data to an output file
    std::string outputFileName = argv[1];
    auto pos = outputFileName.find(".h5");
    if (pos != outputFileName.size()) outputFileName.erase(pos, 3);
    outputFileName += ".Gutzwiller.h5";
    save(outputFileName, mpsApplied);
  }
  return 0;
};