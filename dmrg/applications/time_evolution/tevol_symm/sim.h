/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
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


#include "tevol_sim.h"
#include "tevol_nn_sim.h"
#include "tevol_mpo_sim.h"

template <>
void run_tevol<matrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        tevol_sim<matrix, grp, nearest_neighbors_evolver<matrix, grp> > sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        tevol_sim<matrix, grp, mpo_evolver<matrix, grp> > sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"].str()+")");
    }
}

template <>
void run_tevol<cmatrix, grp>(DmrgParameters & parms, ModelParameters & model)
{
    if (parms["te_type"] == "nn")
    {
        tevol_sim<cmatrix, grp, nearest_neighbors_evolver<cmatrix, grp> > sim(parms, model);
        sim.run();
    } else if (parms["te_type"] == "mpo") {
        tevol_sim<cmatrix, grp, mpo_evolver<cmatrix, grp> > sim(parms, model);
        sim.run();
    } else {
        throw std::runtime_error("Don't know this time evolution. ("+parms["te_type"].str()+")");
    }
}


