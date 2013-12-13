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

#include "dmrg/models/coded/models_u1.hpp"
//#include "dmrg/models/coded/models_bela.hpp"

template<class Matrix>
struct coded_model_factory<Matrix, U1> {
    static boost::shared_ptr<model_impl<Matrix, U1> > parse
    (Lattice const& lattice, BaseParameters & parms)
    {
        typedef boost::shared_ptr<model_impl<Matrix, U1> > impl_ptr;
        if (parms["MODEL"] == std::string("heisenberg"))
            return impl_ptr( new Heisenberg<Matrix>(lattice, parms["Jxy"], parms["Jz"]) );
        else if (parms["MODEL"] == std::string("HCB"))
            return impl_ptr( new HCB<Matrix>(lattice) );
        else if (parms["MODEL"] == std::string("boson Hubbard"))
            return impl_ptr( new BoseHubbard<Matrix>(lattice, parms) );
//        else if (parms["MODEL"] == std::string("fermion Hubbard"))
//            return impl_ptr( new FermiHubbardU1<Matrix>(lattice, parms) );
        else if (parms["MODEL"] == std::string("FreeFermions"))
            return impl_ptr( new FreeFermions<Matrix>(lattice, parms["t"]) );
//        else if (parms["MODEL"] == std::string("bela_chiral"))
//            return impl_ptr( new Chiral<Matrix>(lattice, parms) );
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
