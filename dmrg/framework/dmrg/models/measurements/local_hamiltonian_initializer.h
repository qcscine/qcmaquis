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
    LocalHamiltonianInitializer(std::vector<MPS<Matrix, SymmGroup> >& mps_, // ugly, this should be only MPS
                                MPO<Matrix, SymmGroup> const& mpo_,
                                BaseParameters & parms_,
                                boost::function<bool ()> stop_callback_,
                                int site_) // v---- warning, here we're preparing a copy of our mps that could possibly be avoided
    : base(mps_, mpo_, parms_, stop_callback_, site_, false)
     {
        assert(mps_.size() == 1); // we must pass a vector of size 1 to optimizer_base
     }

    BoundaryVector const & left() {
        // make sure we have only 1 boundary
        assert(this->left_sa_.size() == 1);
        return this->left_sa_[0];
    }
    BoundaryVector const & right() {
        // make sure we have only 1 boundary
        assert(this->right_sa_.size() == 1);
        return this->right_sa_[0];
    }

    // Fetches left and right boundaries from disk for a given site and one/twosite tensor
    // Only intended for left->right sweeps
    void fetch_boundaries(int site, bool twosite = true) {

        std::size_t L = this->mps_vector[0].length();

        assert(this->left_sa_.size() == 1);
        assert(this->right_sa_.size() == 1);

        if (twosite)
            assert(site < L);
        else
            assert(site <= L);

        Storage::prefetch(this->left_sa_[0][site]);
        Storage::fetch(this->left_sa_[0][site]);

        int site2 = twosite ? site+2 : site+1;

        Storage::prefetch(this->right_sa_[0][site2]);
        Storage::fetch(this->right_sa_[0][site2]);


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