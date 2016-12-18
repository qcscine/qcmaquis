/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef GSL_COUPLING_H
#define GSL_COUPLING_H

extern "C" {
    double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc);
    double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf);
    double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji);
}

namespace SU2 {

    inline double mod_coupling(int a, int b, int c,
                        int d, int e, int f,
                        int g, int h, int i)
    {
        double ret = sqrt( (g+1.) * (h+1.) * (c+1.) * (f+1.) ) *
               gsl_sf_coupling_9j(a, b, c,
                                  d, e, f,
                                  g, h, i);
        return ret;
    }

    inline bool triangle(int a, int b, int c)
    {
        return std::abs(a-b) <= c && c <= a+b;
    }

    template <class T>
    inline void set_coupling(int a, int b, int c,
                             int d, int e, int f,
                             int g, int h, int i, T init, T couplings[])
    {
        T prefactor = T(sqrt((i+1.)*(a+1.)/((g+1.)*(c+1.)))) * init;
        if (triangle(a,b,c))
        {
            couplings[0] = prefactor * (T)::SU2::mod_coupling(a, b, c, d, e, f, g, h, i);
            couplings[2] = prefactor * (T)::SU2::mod_coupling(a, b, c, d, e, f, g, 2, i);
        }
        if (triangle(a,2,c))
        {
            couplings[1] = prefactor * (T)::SU2::mod_coupling(a, 2, c, d, e, f, g, h, i);
            couplings[3] = prefactor * (T)::SU2::mod_coupling(a, 2, c, d, e, f, g, 2, i);
        }
    }
}

#endif
