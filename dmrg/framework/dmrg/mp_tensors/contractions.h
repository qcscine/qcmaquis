/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
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

#ifndef CONTRACTIONS_H
#define CONTRACTIONS_H

#ifdef USE_AMBIENT
#include "dmrg/mp_tensors/contractions/impl/ambient.hpp"
#else
//#include "dmrg/mp_tensors/contractions/impl/alps.hpp"
#include "dmrg/mp_tensors/contractions/impl/memsave.hpp"
#endif

#include "dmrg/mp_tensors/contractions/abelian/boundary_times_mps.hpp"

#include "dmrg/mp_tensors/contractions/abelian/apply_op.hpp"
#include "dmrg/mp_tensors/contractions/abelian/move_boundary.hpp"
#include "dmrg/mp_tensors/contractions/abelian/site_hamil.hpp"
#include "dmrg/mp_tensors/contractions/abelian/prediction.hpp"
#include "dmrg/mp_tensors/contractions/abelian/special.hpp"

#include "dmrg/mp_tensors/contractions/abelian_engine.hpp"
#include "dmrg/mp_tensors/contractions/su2_engine.hpp"

#endif
