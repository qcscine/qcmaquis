/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "gsl_coupling.h"

bool WignerWrapper::UseCache = false;

typename WignerWrapper::map_type WignerWrapper::map;

void WignerWrapper::fill_cache(int max)
{
    if (!UseCache) return;

    double dummy;
    for (int i = 0; i <= max; i++)
    for (int j = 0; j <= 2; j++)
    for (int k = 0; k <= max; k++)
    for (int l = 0; l <= 2; l++)
    for (int m = 0; m <= 2; m++)
    for (int n = 0; n <= 2; n++)
    for (int o = 0; o <= max; o++)
    for (int p = 0; p <= 2; p++)
    for (int q = 0; q <= max; q++)
    {
        // maybe not required
        // only include cases l==m==n for l==0
        if ((l == m && m == n) && (l != 0)) continue;

        if (!triangle_9j_fails(i,j,k,l,m,n,o,p,q)) // make sure zero values are not cached
            map[std::make_tuple(i,j,k,l,m,n,o,p,q)] = WignerWrapper::wigner9j_nocache(i,j,k,l,m,n,o,p,q);
    }
}
