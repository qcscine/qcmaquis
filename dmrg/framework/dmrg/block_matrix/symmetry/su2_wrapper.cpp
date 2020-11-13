/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020         Leon Freitag <lefreita@ethz.ch>
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
