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

#include <unordered_map>
#include <mutex>
#include <iostream>
#include "hash_tuple.h"

#include <cmath>

extern "C" {
    double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc);
    double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf);
    double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji);

}
namespace SU2 {
    inline bool triangle(int a, int b, int c)
    {
        return (( a + b + c ) % 2 == 0 ) && std::abs(a-b) <= c && c <= a+b;
    }

}

class WignerWrapper
{
    public:
        // Global variable that enables or disables the cache
        static bool UseCache;

        // Print the cache contents, useful for debugging.
        static void print_map_contents()
        {
            for (auto&& m: map)
                std::cout << std::get<0>(m.first) << " "
                             << std::get<1>(m.first) << " "
                             << std::get<2>(m.first) << " "
                             << std::get<3>(m.first) << " "
                             << std::get<4>(m.first) << " "
                             << std::get<5>(m.first) << " "
                             << std::get<6>(m.first) << " "
                             << std::get<7>(m.first) << " "
                             << std::get<8>(m.first) << " "
                             << m.second << std::endl;
            std::cout << "Number of elements: " << map.size() << std::endl;

        }

        // \brief Fills the Wigner 9j cache with elements up to (max,2,max,2,2,2,max,2,max)
        // Zero elements are not added to the cache
        // \param max maximum index for symbols, calculated as  (max_spin + spin)/2
        //  with max_spin as maximum number of unpaired electrons (see also sim::sim())
        static void fill_cache(int max);

        // woo hoo copy paste!!!
        inline static double gsl_sf_coupling_3j(int two_ja, int two_jb, int two_jc, int two_ma, int two_mb, int two_mc)
        {
            return ::gsl_sf_coupling_3j(two_ja, two_jb, two_jc, two_ma, two_mb, two_mc);
        }

        inline static double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf)
        {

            return ::gsl_sf_coupling_6j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf);
        }
        // \brief Calculate the Wigner 9j symbol, or obtain it from cache
        inline static double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
        {

            double phase = 1.;

            // Consider symmetry properties
            if (triangle_9j_fails(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji)) return 0.0;

            // Reflection along the diagonals does not change the wigner 9j symbol
            if ((two_jb < two_jd) && (two_jc < two_jg) && (two_jf < two_jh))
            {
                std::swap(two_jb, two_jd);
                std::swap(two_jc, two_jg);
                std::swap(two_jf, two_jh);
            }
            else
            if ((two_jb < two_jf) && (two_ja < two_ji) && (two_jd < two_jh))
            {
                std::swap(two_jb, two_jf);
                std::swap(two_ja, two_ji);
                std::swap(two_jd, two_jh);
            }

            // Other Wigner 9j symmetries are ignored because the cache is small enough anyway, and too many symmetry checks here might degrade performance
            // besides, the previous implementations did not seem to work

            return (UseCache) ?
               wigner9j_cache(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji)
                : wigner9j_nocache(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji);
        }


        private:
            typedef std::tuple<int, int, int, int, int, int, int, int, int> gsl_indices;
            typedef std::unordered_map<gsl_indices, double, hash_tuple::hash<gsl_indices> > map_type;

            // The map that stores the values
            static map_type map;

            inline static double wigner9j_cache(int a, int b, int c,
                                    int d, int e, int f,
                                    int g, int h, int i)
            {
                gsl_indices idx=std::make_tuple(a,b,c,d,e,f,g,h,i);
                double ret;
                // Search for the value in cache
                auto map_idx = map.find(idx);
                // If value not found in cache, issue a warning and calculate it
                if (map_idx == map.end())
                {

                    ret = WignerWrapper::wigner9j_nocache(a, b, c,
                                             d, e, f,
                                             g, h, i);

                    // print a warning
                    std::cout << "Warning: Wigner 9j symbol for " << a << "," << b << "," << c <<
                                                              "," << d << "," << e << "," << f <<
                                                              "," << g << "," << h << "," << i << " not found in cache.";

                    // alternatively, add the missing value to the cache, but this is not threadsafe so has been disabled
                    //map[idx] = ret;


                }
                else // use the cached value
                {
                    ret = map_idx->second;
                }
                return ret;
            }

            inline static bool triangle_9j_fails(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
            {
                  return (( !SU2::triangle( two_ja, two_jb, two_jc ) ) ||
                        ( !SU2::triangle( two_jd, two_je, two_jf ) ) ||
                        ( !SU2::triangle( two_jg, two_jh, two_ji ) ) ||
                        ( !SU2::triangle( two_ja, two_jd, two_jg ) ) ||
                        ( !SU2::triangle( two_jb, two_je, two_jh ) ) ||
                        ( !SU2::triangle( two_jc, two_jf, two_ji ) ));

            }
            inline static double wigner9j_nocache(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji)
            {
                return ::gsl_sf_coupling_9j(two_ja, two_jb, two_jc, two_jd, two_je, two_jf, two_jg, two_jh, two_ji);
            }


};

namespace SU2 {

    inline double mod_coupling(int a, int b, int c,
                        int d, int e, int f,
                        int g, int h, int i)
    {
        double ret = sqrt( (g+1.) * (h+1.) * (c+1.) * (f+1.) ) *
               WignerWrapper::gsl_sf_coupling_9j(a, b, c,
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
}

#endif
