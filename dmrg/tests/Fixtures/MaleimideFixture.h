/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2022 Institute for Theoretical Physics, ETH Zurich
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
 *               2022- by Nina Glaser <nglaser@ethz.ch>
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

#ifndef TEST_MALEIMIDE_FIXTURE_H
#define TEST_MALEIMIDE_FIXTURE_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

/** @brief Fixture class for the test of the Watson Hamiltonian of maleimide */
struct MaleimideFixture
{
  // Types definition
  using MaquisIntegralType = maquis::integral_map<double, chem::Hamiltonian::VibrationalCanonical>;
  /** @brief Constructor for the fixture class */
  MaleimideFixture() {
    integralFileOneBody.open("integral_file_test_Watson_Maleimide_OneBody");
    integralFileOneBody << " 3.16758229e+01  1  1  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.16758229e+01 -1 -1  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.54421822e+00  1  1  1  1  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.27214847e-02  1  1  1  1  1  1  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.02901434e-04  1  1  1  1  1  1  1  1  0  0  0  0" << std::endl;
    integralFileOneBody << " 6.87503072e-06  1  1  1  1  1  1  1  1  1  1  0  0" << std::endl;
    integralFileOneBody << "-1.08654847e-07  1  1  1  1  1  1  1  1  1  1  1  1" << std::endl;
    integralFileOneBody << " 7.32647784e+01  2  2  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-7.32647784e+01 -2 -2  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.08258673e+00  2  2  2  2  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-7.56289543e-03  2  2  2  2  2  2  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.50264835e-04  2  2  2  2  2  2  2  2  0  0  0  0" << std::endl;
    integralFileOneBody << "-5.03781540e-06  2  2  2  2  2  2  2  2  2  2  0  0" << std::endl;
    integralFileOneBody << " 7.90742573e-08  2  2  2  2  2  2  2  2  2  2  2  2" << std::endl;
    integralFileOneBody << " 1.04890799e+00  3  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 9.83195240e+01  3  3  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-9.83195240e+01 -3 -3  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.44388799e-01  3  3  3  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.15799869e-01  3  3  3  3  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.33180702e-04  3  3  3  3  3  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.94391618e-04  3  3  3  3  3  3  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.92899693e-08  3  3  3  3  3  3  3  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.69005265e-07  3  3  3  3  3  3  3  3  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.39992878e-09  3  3  3  3  3  3  3  3  3  3  0  0" << std::endl;
    integralFileOneBody << " 1.24582125e+02  4  4  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.24582125e+02 -4 -4  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.54524174e+01  4  4  4  4  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.36498385e+00  4  4  4  4  4  4  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.42579007e-02  4  4  4  4  4  4  4  4  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.10610803e-03  4  4  4  4  4  4  4  4  4  4  0  0" << std::endl;
    integralFileOneBody << " 2.67158541e-05  4  4  4  4  4  4  4  4  4  4  4  4" << std::endl;
    integralFileOneBody << " 1.34568908e+02  5  5  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.34568908e+02 -5 -5  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.55386274e-01  5  5  5  5  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.60347444e-04  5  5  5  5  5  5  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.37134929e-06  5  5  5  5  5  5  5  5  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.12024367e-07  5  5  5  5  5  5  5  5  5  5  0  0" << std::endl;
    integralFileOneBody << "-1.67437069e-09  5  5  5  5  5  5  5  5  5  5  5  5" << std::endl;
    integralFileOneBody << " 1.57549454e+02  6  6  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.57549454e+02 -6 -6  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.97079482e+00  6  6  6  6  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.76246912e-02  6  6  6  6  6  6  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.50808429e-04  6  6  6  6  6  6  6  6  0  0  0  0" << std::endl;
    integralFileOneBody << "-6.07447689e-06  6  6  6  6  6  6  6  6  6  6  0  0" << std::endl;
    integralFileOneBody << " 8.05542890e-08  6  6  6  6  6  6  6  6  6  6  6  6" << std::endl;
    integralFileOneBody << " 6.05542440e-01  7  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.59809793e+02  7  7  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.59809793e+02 -7 -7  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.01240047e+00  7  7  7  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.37998118e-02  7  7  7  7  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.55630068e-05  7  7  7  7  7  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.28676299e-05  7  7  7  7  7  7  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.35749440e-06  7  7  7  7  7  7  7  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.20993831e-06  7  7  7  7  7  7  7  7  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.79155795e-08  7  7  7  7  7  7  7  7  7  0  0  0" << std::endl;
    integralFileOneBody << " 1.59963018e-07  7  7  7  7  7  7  7  7  7  7  0  0" << std::endl;
    integralFileOneBody << "-2.86768012e-10  7  7  7  7  7  7  7  7  7  7  7  0" << std::endl;
    integralFileOneBody << "-2.25040325e-09  7  7  7  7  7  7  7  7  7  7  7  7" << std::endl;
    integralFileOneBody << " 1.68916764e+02  8  8  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.68916764e+02 -8 -8  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 5.82129849e-02  8  8  8  8  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.23422939e-04  8  8  8  8  8  8  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 7.56145654e-07  8  8  8  8  8  8  8  8  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.25029381e-08  8  8  8  8  8  8  8  8  8  8  0  0" << std::endl;
    integralFileOneBody << " 6.44667508e-10  8  8  8  8  8  8  8  8  8  8  8  8" << std::endl;
    integralFileOneBody << " 1.94184345e+02  9  9  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.94184345e+02 -9 -9  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 7.38993626e-02  9  9  9  9  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.52236224e-03  9  9  9  9  9  9  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.10894646e-04  9  9  9  9  9  9  9  9  0  0  0  0" << std::endl;
    integralFileOneBody << "-7.21463486e-06  9  9  9  9  9  9  9  9  9  9  0  0" << std::endl;
    integralFileOneBody << " 9.45752150e-08  9  9  9  9  9  9  9  9  9  9  9  9" << std::endl;
    integralFileOneBody << " 2.09253105e+02  10  10  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.09253105e+02 -10 -10  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.88335711e+00  10  10  10  10  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.76545820e-02  10  10  10  10  10  10  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.95950286e-04  10  10  10  10  10  10  10  10  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.95858267e-05  10  10  10  10  10  10  10  10  10  10  0  0" << std::endl;
    integralFileOneBody << " 4.22219971e-07  10  10  10  10  10  10  10  10  10  10  10  10" << std::endl;
    integralFileOneBody << " 1.19353309e-01  11  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.27618685e+02  11  11  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.27618685e+02 -11 -11  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.78018793e+00  11  11  11  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.28246780e-01  11  11  11  11  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.70625801e-03  11  11  11  11  11  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 7.32642995e-05  11  11  11  11  11  11  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.94486113e-06  11  11  11  11  11  11  11  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-7.78857466e-10  11  11  11  11  11  11  11  11  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.42293770e-07  11  11  11  11  11  11  11  11  11  0  0  0" << std::endl;
    integralFileOneBody << "-1.21453210e-08  11  11  11  11  11  11  11  11  11  11  0  0" << std::endl;
    integralFileOneBody << " 3.47029232e-09  11  11  11  11  11  11  11  11  11  11  11  0" << std::endl;
    integralFileOneBody << " 4.23237194e-10  11  11  11  11  11  11  11  11  11  11  11  11" << std::endl;
    integralFileOneBody << " 2.30513204e+02  12  12  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.30513204e+02 -12 -12  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.27591067e-01  12  12  12  12  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-9.17475164e-04  12  12  12  12  12  12  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 7.21116110e-06  12  12  12  12  12  12  12  12  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.29400879e-07  12  12  12  12  12  12  12  12  12  12  0  0" << std::endl;
    integralFileOneBody << " 1.24529042e-09  12  12  12  12  12  12  12  12  12  12  12  12" << std::endl;
    integralFileOneBody << " 2.40978301e+02  13  13  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.40978301e+02 -13 -13  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.75162373e+00  13  13  13  13  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.99738879e-02  13  13  13  13  13  13  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.02496262e-03  13  13  13  13  13  13  13  13  0  0  0  0" << std::endl;
    integralFileOneBody << "-9.01165984e-05  13  13  13  13  13  13  13  13  13  13  0  0" << std::endl;
    integralFileOneBody << " 1.20625613e-06  13  13  13  13  13  13  13  13  13  13  13  13" << std::endl;
    integralFileOneBody << "-1.53831111e+00  14  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.68942645e+02  14  14  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.68942645e+02 -14 -14  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.17853165e+00  14  14  14  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.73288644e+00  14  14  14  14  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.19579769e-02  14  14  14  14  14  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.65012679e-02  14  14  14  14  14  14  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.02318851e-04  14  14  14  14  14  14  14  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 9.48561644e-05  14  14  14  14  14  14  14  14  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.58077261e-07  14  14  14  14  14  14  14  14  14  0  0  0" << std::endl;
    integralFileOneBody << "-2.61085518e-07  14  14  14  14  14  14  14  14  14  14  0  0" << std::endl;
    integralFileOneBody << "-9.20028677e-09  14  14  14  14  14  14  14  14  14  14  14  0" << std::endl;
    integralFileOneBody << "-6.83475836e-10  14  14  14  14  14  14  14  14  14  14  14  14" << std::endl;
    integralFileOneBody << " 2.88856434e+02  15  15  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.88856434e+02 -15 -15  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.60050669e-01  15  15  15  15  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-7.22110426e-04  15  15  15  15  15  15  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.59337302e-06  15  15  15  15  15  15  15  15  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.96353933e-07  15  15  15  15  15  15  15  15  15  15  0  0" << std::endl;
    integralFileOneBody << " 2.33061941e-09  15  15  15  15  15  15  15  15  15  15  15  15" << std::endl;
    integralFileOneBody << " 3.32692453e+02  16  16  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.32692453e+02 -16 -16  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.30154893e-01  16  16  16  16  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-2.80683467e-03  16  16  16  16  16  16  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.47398660e-05  16  16  16  16  16  16  16  16  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.11547756e-07  16  16  16  16  16  16  16  16  16  16  0  0" << std::endl;
    integralFileOneBody << " 6.64381376e-10  16  16  16  16  16  16  16  16  16  16  16  16" << std::endl;
    integralFileOneBody << " 3.37429340e+02  17  17  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.37429340e+02 -17 -17  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.17001398e+00  17  17  17  17  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.78182725e-02  17  17  17  17  17  17  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.75955813e-04  17  17  17  17  17  17  17  17  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.93142458e-06  17  17  17  17  17  17  17  17  17  17  0  0" << std::endl;
    integralFileOneBody << " 1.69327600e-08  17  17  17  17  17  17  17  17  17  17  17  17" << std::endl;
    integralFileOneBody << " 3.50821238e-01  18  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.40425861e+02  18  18  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.40425861e+02 -18 -18  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.90191694e+00  18  18  18  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.66368953e-03  18  18  18  18  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.86036997e-03  18  18  18  18  18  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.67743403e-04  18  18  18  18  18  18  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.69367040e-07  18  18  18  18  18  18  18  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.61601098e-06  18  18  18  18  18  18  18  18  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.20384240e-08  18  18  18  18  18  18  18  18  18  0  0  0" << std::endl;
    integralFileOneBody << "-6.91309813e-08  18  18  18  18  18  18  18  18  18  18  0  0" << std::endl;
    integralFileOneBody << " 1.64322543e-09  18  18  18  18  18  18  18  18  18  18  18  0" << std::endl;
    integralFileOneBody << " 1.07895954e-09  18  18  18  18  18  18  18  18  18  18  18  18" << std::endl;
    integralFileOneBody << "-2.79273844e-01  19  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.04469600e+02  19  19  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.04469600e+02 -19 -19  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.47869229e+01  19  19  19  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.20580495e+00  19  19  19  19  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.67307777e-02  19  19  19  19  19  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.51811460e-03  19  19  19  19  19  19  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.67405642e-05  19  19  19  19  19  19  19  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.58884557e-06  19  19  19  19  19  19  19  19  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.51566353e-08  19  19  19  19  19  19  19  19  19  0  0  0" << std::endl;
    integralFileOneBody << "-4.05745367e-09  19  19  19  19  19  19  19  19  19  19  0  0" << std::endl;
    integralFileOneBody << "-8.15571255e-10  19  19  19  19  19  19  19  19  19  19  19  0" << std::endl;
    integralFileOneBody << " 4.50612833e+02  20  20  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.50612833e+02 -20 -20  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.19801174e-01  20  20  20  20  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.66376722e-04  20  20  20  20  20  20  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.29003628e-05  20  20  20  20  20  20  20  20  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.99161789e-06  20  20  20  20  20  20  20  20  20  20  0  0" << std::endl;
    integralFileOneBody << " 4.05785099e-08  20  20  20  20  20  20  20  20  20  20  20  20" << std::endl;
    integralFileOneBody << "-6.14602515e-02  21  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.59825663e+02  21  21  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.59825663e+02 -21 -21  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.18841856e+01  21  21  21  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.64883340e-01  21  21  21  21  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 2.62814679e-02  21  21  21  21  21  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 6.78162378e-04  21  21  21  21  21  21  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.80997855e-05  21  21  21  21  21  21  21  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 3.35648222e-07  21  21  21  21  21  21  21  21  0  0  0  0" << std::endl;
    integralFileOneBody << "-1.32648603e-07  21  21  21  21  21  21  21  21  21  0  0  0" << std::endl;
    integralFileOneBody << "-1.59050576e-08  21  21  21  21  21  21  21  21  21  21  0  0" << std::endl;
    integralFileOneBody << " 2.82400172e-09  21  21  21  21  21  21  21  21  21  21  21  0" << std::endl;
    integralFileOneBody << " 3.74621008e-10  21  21  21  21  21  21  21  21  21  21  21  21" << std::endl;
    integralFileOneBody << " 8.09908987e+02  22  22  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.09908987e+02 -22 -22  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 5.60811556e+00  22  22  22  22  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.56026384e-02  22  22  22  22  22  22  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 6.66358446e-05  22  22  22  22  22  22  22  22  0  0  0  0" << std::endl;
    integralFileOneBody << "-4.03996375e-03  23  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.14839990e+02  23  23  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.14839990e+02 -23 -23  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.11436328e+01  23  23  23  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 5.50948413e+00  23  23  23  23  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.11763528e-01  23  23  23  23  23  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.53014547e-02  23  23  23  23  23  23  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-8.56837243e-04  23  23  23  23  23  23  23  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 7.37258952e-05  23  23  23  23  23  23  23  23  0  0  0  0" << std::endl;
    integralFileOneBody << "-3.16883574e-06  23  23  23  23  23  23  23  23  23  0  0  0" << std::endl;
    integralFileOneBody << "-3.42133088e-07  23  23  23  23  23  23  23  23  23  23  0  0" << std::endl;
    integralFileOneBody << " 2.67639062e-08  23  23  23  23  23  23  23  23  23  23  23  0" << std::endl;
    integralFileOneBody << " 6.86217417e-10  23  23  23  23  23  23  23  23  23  23  23  23" << std::endl;
    integralFileOneBody << " 2.81543538e-01  24  0  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 9.15244681e+02  24  24  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << "-9.15244681e+02 -24 -24  0  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.35512484e+02  24  24  24  0  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.35964669e+01  24  24  24  24  0  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 1.13306417e+00  24  24  24  24  24  0  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 8.23025023e-02  24  24  24  24  24  24  0  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 6.14190147e-03  24  24  24  24  24  24  24  0  0  0  0  0" << std::endl;
    integralFileOneBody << " 5.69417505e-04  24  24  24  24  24  24  24  24  0  0  0  0" << std::endl;
    integralFileOneBody << " 4.46980482e-05  24  24  24  24  24  24  24  24  24  0  0  0" << std::endl;
    integralFileOneBody << " 1.41467672e-06  24  24  24  24  24  24  24  24  24  24  0  0" << std::endl;
    integralFileOneBody << "-4.63710151e-08  24  24  24  24  24  24  24  24  24  24  24  0" << std::endl;
    integralFileOneBody << "-2.78835504e-09  24  24  24  24  24  24  24  24  24  24  24  24" << std::endl;
    integralFileOneBody.close();
    // Sets the parameters
    parametersMaleimideOneBody.set("L", 24);
    parametersMaleimideOneBody.set("symmetry", "none");
    parametersMaleimideOneBody.set("LATTICE", "watson lattice");
    parametersMaleimideOneBody.set("MODEL", "watson");
    parametersMaleimideOneBody.set("Nmax", "2,2,3,2,3,2,3,3,3,2,3,3,2,2,3,3,2,3,3,2,3,2,3,4");
    parametersMaleimideOneBody.set("integral_file", "integral_file_test_Watson_Maleimide_OneBody");
    parametersMaleimideOneBody.set("watson_coordinate_type", "cartesian");
    parametersMaleimideOneBody.set("watson_max_coupling_input", 12);
    parametersMaleimideOneBody.set("watson_max_coupling", 1);
  }

  /** @brief Class destructor (removes tmp files) */
  ~MaleimideFixture() {
    std::remove("integral_file_test_Watson_Maleimide_OneBody");
  }

  // Class members
  DmrgParameters parametersMaleimideOneBody;
  std::ofstream integralFileOneBody;
};

#endif // TEST_MALEIMIDE_FIXTURE_H