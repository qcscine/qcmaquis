/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef CONTRACTIONS_SU2_TEST_COUPLING_H
#define CONTRACTIONS_SU2_TEST_COUPLING_H

#include "dmrg/mp_tensors/contractions/non-abelian/gsl_coupling.h"

namespace SU2 {

    double destroy(int jR, int jRt, int local_spin) {
        // jR = bra right spin, jRt = ket right spin

        double ret = 1.;
        if (local_spin == 0) {
            ret = std::sqrt( (jR + 1.) / (jRt + 1.));
            int phase = (jRt - jR + 1)/2;
            for (int p_= 0; p_ < phase; ++p_) ret *= -1.;
        }
        return ret;
    }

    double F0_create(int jR, int jLt, int local_spin, int spin_out) {
        // jR = bra..
        double ret;
        if(local_spin == 0) {
            ret = 1./sqrt(2.);
        }
        else {
            ret = std::sqrt( (jLt + 1.) / (jR + 1.));
            ret *= 1./sqrt(2.);
            int phase = (jR - jLt + 1)/2;
            for (int p_= 0; p_ < phase; ++p_) ret *= -1.;
        }
        return ret;
    }

    double S0c(int S1, int SLU, int SLD, int SD, int NU, int NLU) {
        //std::swap(SLU,SLD);
        //std::swap(S1,SD);
        if (NU-NLU == 0) {
            int fase = ((((S1 - SLD + 1)/2)%2)!=0)?-1:1;
            return fase * sqrt(0.5 * (SLD+1.0) / (S1+1.0) );
        }
        else {
            if (NU-NLU != 1) return 9999;
            return -sqrt(0.5);
        }
    }

    double S1c(int S1, int SLU, int SLD, int SD, int NU, int NLU) {
        //std::swap(SLU,SLD);
        //std::swap(S1,SD);
        if (NU-NLU == 0) {
            int fase = ((((S1 + SD + 2)/2)%2)!=0)?-1:1;
            return fase * sqrt(3.0*(SLD+1)) * gsl_sf_coupling_6j(1,1,2,S1,SD,SLD);
        }
        else {
            if (NU-NLU != 1) return 9999;
            int fase = ((((SLU + SD + 1)/2)%2)!=0)?-1:1;
            return fase * sqrt(3.0*(S1+1)) * gsl_sf_coupling_6j(1,1,2,S1,SD,SLU);
        }
    }
}

#endif
