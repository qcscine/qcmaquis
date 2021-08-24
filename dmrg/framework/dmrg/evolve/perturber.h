/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef PERTURBER_H
#define PERTURBER_H

#include "dmrg/mp_tensors/mpstensor.h"
#include "dmrg/mp_tensors/boundary.h"
#include "dmrg/mp_tensors/mpo.h"

/**
 * @brief Class managing the addition of the noise.
 * 
 * For now, this class implements only the noise model presented in the paper
 * "A Strictly Single-Site DMRG Algorithm with Subspace Expansion"
 *  Phys. Rev. B 91, 155115 (2015)
 * In the future, it will support also density matrix-based perturbations.
 */

template<class Matrix, class SymmGroup>
class Perturber
{
    // Types definition
    using boundary_type = std::vector< Boundary<Matrix, SymmGroup> >;
    using alpha_type = double;
    using MPSTensor_type = MPSTensor<Matrix, SymmGroup>;
    using MPO_type = MPO<Matrix, SymmGroup>;
    using site_type = std::size_t;
public:
    // Constructor
    template<class BaseParameters>
    Perturber(boundary_type const& left, boundary_type const& right, MPO_type const& MPO, BaseParameters& input_parms) ;
    // Setter
    void activate_perturb_dm() { perturb_dm_ = true ; };
    void activate_perturb_mps() { perturb_MPS_ = true ; };
    void update_alpha(alpha_type alpha_new);
    // Getter
    bool is_perturb_mps() const { return this->perturb_MPS_ ; };
    bool is_perturb_dm() const { return this->perturb_dm_; };
    auto getAlpha() const { return this->alpha_; };
    // Relevant methods
    void perturb_mps(MPSTensor_type& MPS, site_type site, bool is_left) const;
private:
    // -- Private attributes --
    const boundary_type& left_bound_;
    const boundary_type& right_bound_;
    const MPO_type& MPO_;
    bool perturb_dm_, perturb_MPS_;
    alpha_type alpha_, thresh_scal_;
};

#include "perturber.cpp"

#endif