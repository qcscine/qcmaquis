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
#ifndef LOCAL_HAMILTONIAN_INITIALIZER_H
#define LOCAL_HAMILTONIAN_INITIALIZER_H

#include "dmrg/optimize/optimize.h"

// LocalHamiltonianInitializer
// Constructs boundaries used in the calculation of the local Hamiltonian matrix elements
template<class Matrix, class SymmGroup, class Storage>
class LocalHamiltonianInitializer : public optimizer_base<Matrix, SymmGroup, Storage>
{
    typedef optimizer_base<Matrix, SymmGroup, Storage> base;
    typedef std::vector<Boundary<typename storage::constrained<Matrix>::type, SymmGroup> > BoundaryVector;

    public:
    LocalHamiltonianInitializer(MPS<Matrix, SymmGroup> & mps_,
                                MPO<Matrix, SymmGroup> const& mpo_,
                                BaseParameters & parms_,
                                boost::function<bool ()> stop_callback_,
                                int site_)
    : base(mps_, mpo_, parms_, stop_callback_, site_, false) {}

    BoundaryVector const & left() { return this->left_; }
    BoundaryVector const & right() { return this->right_; }

    // Fetches left and right boundaries from disk for a given site and one/twosite tensor
    // Only intended for left->right sweeps, i.e. the boundaries for site L-1 are fetched at site L-2
    void fetch_boundaries(int site, bool twosite = true) {

        std::size_t L = this->mps.length();

        assert(site < L);

        Storage::prefetch(this->left_[site]);
        Storage::fetch(this->left_[site]);

        if (twosite)
        {
            int site2 = site < L-1 ? site+2 : site+1;

            Storage::prefetch(this->right_[site2]);
            Storage::fetch(this->right_[site2]);
        }
        else
        {
            Storage::prefetch(this->right_[site+1]);
            Storage::fetch(this->right_[site+1]);
        }

    }

    // Removes (unneeded) boundaries from memory
    // TODO: Figure out how to do it properly!
    void drop_boundaries(int site, bool twosite = true) {
        throw std::runtime_error("LocalHamiltonianInitializer::drop_boundaries() not implemented yet");
    }

    virtual void sweep(int sweep, OptimizeDirection d = Both) {
        throw std::runtime_error("LocalHamiltonianInitializer::sweep() should not be called");
    }
};
#endif