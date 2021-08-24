/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2018-2020 by Alberto Baiardi <abaiardi@ethz.ch>
 *
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 2 or (at your option) any later
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

#include <chrono>
#include "perturber.h"
#include "dmrg/mp_tensors/contractions.h"

// -- CONSTRUCTOR --
template<class Matrix, class SymmGroup>
template<class BaseParameters>
Perturber<Matrix, SymmGroup>::Perturber(Perturber::boundary_type const& left, Perturber::boundary_type const& right,
                                        Perturber::MPO_type const& MPO, BaseParameters& input_pars)
    : left_bound_(left), right_bound_(right), MPO_(MPO), alpha_(0.), perturb_dm_(false),
      perturb_MPS_(false), thresh_scal_(1.0E-3)
{
    if (input_pars.is_set("perturbation_type")) {
        if (input_pars["perturbation_type"] == "mps")
            activate_perturb_mps();
        else if (input_pars["perturbation_type"] == "density_matrix")
            activate_perturb_dm();
    }
}

template<class Matrix, class SymmGroup>
void Perturber<Matrix, SymmGroup>::update_alpha(Perturber::alpha_type alpha_new)
{
    alpha_ = alpha_new ;
    if (alpha_ < 1.0E-30)
        perturb_MPS_ = false ;
}

template<class Matrix, class SymmGroup>
void Perturber<Matrix, SymmGroup>::perturb_mps(MPSTensor_type& MPS, site_type site, bool is_left) const
{
    // Builds the perturber (which is a boundary)
    typename Boundary<Matrix, SymmGroup>::Boundary perturber;
    // Generates the correct pairing
    if (perturb_MPS_) {
        if (is_left) {
            MPS.make_left_paired();
            perturber = contraction::Engine<Matrix, Matrix, SymmGroup>::generate_left_mpo_basis(MPS, MPS, left_bound_[site], MPO_[site]);
            perturber *= alpha_;
            // Merges all the contributions of the perturber
            auto ToAddBM = block_matrix<Matrix, SymmGroup>();
            for (std::size_t b = 0; b < perturber.aux_dim(); ++b) {
                for (std::size_t idx = 0; idx < perturber[b].n_blocks(); idx++)
                    if (perturber[b].get_left_charge(idx) == perturber[b].get_right_charge(idx) && perturber[b].norm() > 1.0E-10)
                            ToAddBM.add_block_to_column(perturber[b], perturber[b].get_left_charge(idx),
                                                        perturber[b].get_right_charge(idx));
            }
            if (ToAddBM.n_blocks() != 0)
                MPS.add_block_to_column(ToAddBM);
        } else {
            MPS.make_right_paired();
            perturber = contraction::Engine<Matrix, Matrix, SymmGroup>::generate_right_mpo_basis(MPS, MPS, right_bound_[site+1], MPO_[site]);
            perturber *= alpha_;
            auto ToAddBM = block_matrix<Matrix, SymmGroup>();
            for (int b = 0; b < perturber.aux_dim(); ++b)
                for (int idx = 0; idx < perturber[b].n_blocks(); idx++)
                    if (perturber[b].get_left_charge(idx) == perturber[b].get_right_charge(idx) && perturber[b].norm() > 1.0E-10)
                        ToAddBM.add_block_to_row(perturber[b], perturber[b].get_left_charge(idx),
                                                 perturber[b].get_right_charge(idx));
            if (ToAddBM.n_blocks() != 0)
                MPS.add_block_to_row(ToAddBM);
        }
    }
}

