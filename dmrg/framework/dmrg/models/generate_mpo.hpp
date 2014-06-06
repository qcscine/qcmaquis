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

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "dmrg/models/generate_mpo/mpo_maker.hpp"
#include "dmrg/models/generate_mpo/tagged_mpo_maker_optim.hpp"
#include "dmrg/models/generate_mpo/corr_maker.hpp"
#include "dmrg/models/generate_mpo/bg_corr_maker.hpp"

#include "dmrg/models/chem/pg_symm_converter.h"

#include "dmrg/models/model.h"


template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> make_mpo(Lattice const& lat, Model<Matrix, SymmGroup> const& model, BaseParameters& parms)
{
    //    typename Model<Matrix, SymmGroup>::terms_type const& terms = model.hamiltonian_terms();
    //    generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpom(lat.size(), model.identity_matrix_tag(), model.operators_table());
    //    for (std::size_t i = 0; i < terms.size(); ++i)
    //        mpom.add_term(terms[i], model.filling_matrix_tag());
    
    generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpom(lat, model);
    MPO<Matrix, SymmGroup> mpo = mpom.create_mpo();
    
    PGSymmetryConverter<Matrix, SymmGroup> symm_conv(lat);
    symm_conv.convert_tags_to_symm_tags(mpo);
    
    typedef typename MPOTensor<Matrix, SymmGroup>::col_proxy col_proxy;
    typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
    for (int p = 0; p < lat.size(); ++p) {
        for (int mpocol = 0; mpocol < mpo[p].col_dim(); ++mpocol) {
            col_proxy col_b2 = mpo[p].column(mpocol);
            for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
                index_type b1 = col_it.index();
                MPOTensor_detail::term_descriptor<Matrix, SymmGroup> access = mpo[p].at(b1,mpocol);
                block_matrix<Matrix, SymmGroup> const & W = access.op;
                
                maquis::cout << "site " << p << "(" << b1 << "," << mpocol << ")\n" << W << std::endl;
            }
        }
    }
    
    return mpo;
}

#endif
