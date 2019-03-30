/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#include <unordered_map>
#include <mutex>
#include "hash_tuple.h"

extern "C" {
    double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc);
    double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf);
    double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji);
}

namespace SU2 {

    inline bool triangle(int a, int b, int c)
    {
        return std::abs(a-b) <= c && c <= a+b;
    }
    // Cache for Wigner-9j symbols
    // Profilers show that over 60% of the time is spent in gsl_sf_coupling_9j
    // therefore we need a class that caches the calls to this function
    // uses std::unordered map as a cache
    // Currently we don't exploit the permutation symmetry of the Wigner 9j symbols yet
    // but this should be done in the future
    // The cache should be filled with fill_cache() as early as possible when the parameters
    // (spin, L, nelec) are already known. Unfortunately it's not possible to know them
    // in the constructor already, so we have to use such a kludge.
    class GSLCouplingCache
    {
        private:
            typedef std::tuple<int, int, int, int, int, int, int, int, int> gsl_indices;
            // derive gsl_indices from tuple to override operator== to incorporate symmetry

            // The map that stores the values
            std::unordered_map<gsl_indices, double, hash_tuple::hash<gsl_indices> > map;

        public:
            // Fill the Wigner 9j cache
            // Ideally this should happen in some constructor, but as a static variable
            // the cache has no knowledge about the symmetries and the spin
            // so we need to run fill_cache manually as early as possible in the simulation
            // TODO: protect it better
            void fill_cache(int nelec, int L, int spin)
            {
                int max_idx = std::max((L - std::abs(nelec - L) + spin)/2, 1);
                int max_idx_i = std::max(max_idx, 2);
                std::mutex m_;

                for (int a = 0; a <= max_idx; a++)
                for (int b = 0; b <= 2; b++)
                for (int c = 0; c <= max_idx; c++)
                for (int d = 0; d <= 2; d++)
                for (int e = 0; e <= 2; e++)
                for (int f = 0; f <= 2; f++)
                for (int g = 0; g <= max_idx; g++)
                for (int h = 0; h <= 2; h++)
                for (int i = 0; i <= max_idx_i; i++)
                {
                    gsl_indices idx = std::make_tuple(a,b,c,d,e,f,g,h,i);

                    std::lock_guard<std::mutex> lk(m_);
                    map[idx] = sqrt( (g+1.) * (h+1.) * (c+1.) * (f+1.) ) *
                               gsl_sf_coupling_9j(a,b,c,d,e,f,g,h,i);
                }
            }
            inline double mod_coupling(int a, int b, int c,
                                    int d, int e, int f,
                                    int g, int h, int i)
            {
                gsl_indices idx = std::make_tuple(a,b,c,d,e,f,g,h,i);
                return map[idx];
                // if the map hasn't been initialized this will fail!
                // But we avoid exception handling for the performance
            }

            template <class T>
            inline void set_coupling(int a, int b, int c,
                                     int d, int e, int f,
                                     int g, int h, int i, T init, T couplings[])
            {
                T prefactor = T(sqrt((i+1.)*(a+1.)/((g+1.)*(c+1.)))) * init;
                if (triangle(a,b,c))
                {
                    couplings[0] = prefactor * (T)mod_coupling(a, b, c, d, e, f, g, h, i);
                    couplings[2] = prefactor * (T)mod_coupling(a, b, c, d, e, f, g, 2, i);
                }
                if (triangle(a,2,c))
                {
                    couplings[1] = prefactor * (T)mod_coupling(a, 2, c, d, e, f, g, h, i);
                    couplings[3] = prefactor * (T)mod_coupling(a, 2, c, d, e, f, g, 2, i);
                }
            }

    };

    static GSLCouplingCache gsl_coupling_cache;

    inline double mod_coupling(int a, int b, int c,
                        int d, int e, int f,
                        int g, int h, int i)
    {
        return gsl_coupling_cache.mod_coupling(a,b,c,d,e,f,g,h,i);
    }

    template <class T>
    inline void set_coupling(int a, int b, int c,
                             int d, int e, int f,
                             int g, int h, int i, T init, T couplings[])
    {
        gsl_coupling_cache.set_coupling<T>(a,b,c,d,e,f,g,h,i,init,couplings);
    }

    // No cache
/*
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

*/
}

#endif
