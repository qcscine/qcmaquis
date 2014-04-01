/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_IDENTITY_MPS_H
#define MAQUIS_DMRG_IDENTITY_MPS_H

#include "dmrg/block_matrix/grouped_symmetry.h"

template <class Matrix, class InSymm>
MPS<Matrix, typename grouped_symmetry<InSymm>::type> identity_dm_mps(std::size_t L, Index<InSymm> const& phys_psi,
                                                                     std::vector<Index<typename grouped_symmetry<InSymm>::type> > const& allowed)
{
    MPOTensor<Matrix, InSymm> t(1,1);
    t.set(0, 0, identity_matrix<Matrix>(phys_psi));

    MPO<Matrix, InSymm> mpo(L, t);
    return mpo_to_smps_group(mpo, phys_psi, allowed);
}

#endif
