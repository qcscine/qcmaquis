/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef GENERATE_MPO_H
#define GENERATE_MPO_H

#include "dmrg/models/generate_mpo/mpo_maker.hpp"
#include "dmrg/models/generate_mpo/tagged_mpo_maker_optim.hpp"
#include "dmrg/models/generate_mpo/corr_maker.hpp"

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
    
    PGSymmetryConverter<Matrix, SymmGroup> symm_conv( parse_symm<SymmGroup>(lat.size(), parms) );
    symm_conv.convert_tags_to_symm_tags(mpo);
    
    return mpo;
}

#endif
