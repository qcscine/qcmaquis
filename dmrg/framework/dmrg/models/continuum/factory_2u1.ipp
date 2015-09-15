/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#include "dmrg/models/continuum/models_2u1.hpp"
//#include "dmrg/models/continuum/super_models_2u1.hpp"

template<class Matrix>
struct cont_model_factory<Matrix, TwoU1> {
    static boost::shared_ptr<model_impl<Matrix, TwoU1> > parse
    (Lattice const & lattice, BaseParameters & model)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TwoU1> > impl_ptr;
        std::string model_str = model.is_set("model") ? "model" : "MODEL";
        if (model[model_str] == std::string("fermi_optical_lattice"))
            return impl_ptr( new FermiOpticalLattice<Matrix>(lattice, model) );
//        else if (model[model_str] == std::string("optical_lattice_cons_dm"))
//            return impl_ptr( new DMOpticalLatticeTwoU1<Matrix>(lattice, model) );
//        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
//        }
    }
};
