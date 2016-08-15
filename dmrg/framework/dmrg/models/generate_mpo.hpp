/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "dmrg/models/generate_mpo/mpo_maker.hpp"
#include "dmrg/models/generate_mpo/tagged_mpo_maker_optim.hpp"
#include "dmrg/models/generate_mpo/corr_maker.hpp"
#include "dmrg/models/generate_mpo/1D_mpo_maker.hpp"


template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_mpo(Lattice const& lat, Model<Matrix, SymmGroup> & model)
{
    model.create_terms();
    generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpom(lat, model);
    MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();

    return mpo;
}

#endif
