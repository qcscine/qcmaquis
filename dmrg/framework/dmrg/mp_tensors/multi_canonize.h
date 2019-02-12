/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef MULTI_CANONIZE_H
#define MULTI_CANONIZE_H
/*
    Here we implement a simultaneous canonicalisation of MPS that have been
    optimised simultaneously
*/
#include "dmrg/mp_tensors/mpstensor.h"
#include <functional>

template<class Matrix, class SymmGroup> class MPS;
// multi_normalize_left()
// Split a vector of MPSTensors into a single left-normalised MPSTensor and a vector of MPSTensor
// which can be premultiplied by MPSTensors at the next site to obtain a mixed-canonical form of
// several MPS with the same basis
template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > multi_normalize_left(std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > >& mps_vec, DecompMethod method);

template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > multi_normalize_right(std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > >& mps_vec, DecompMethod method);

// Moves normalisation of an MPS vector to the right, just as the single-mps version in mps.hpp
template<class Matrix, class SymmGroup>
void multi_move_normalization_l2r(std::vector<MPS<Matrix, SymmGroup> > & vec, size_t p1, size_t p2, DecompMethod method);

template<class Matrix, class SymmGroup>
void multi_move_normalization_r2l(std::vector<MPS<Matrix, SymmGroup> > & vec, size_t p1, size_t p2, DecompMethod method);

template<class Matrix, class SymmGroup>
void multi_canonize(std::vector<MPS<Matrix, SymmGroup> > & vec, std::size_t center, DecompMethod method = QR);

#include "dmrg/mp_tensors/multi_canonize.hpp"
#endif