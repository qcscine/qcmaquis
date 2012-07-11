/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifdef USE_AMBIENT
// parallel matrix
#include "dmrg/block_matrix/detail/ambient.hpp"
#endif

#include <complex>
#include <vector>

// serial matrix
#include "dmrg/block_matrix/detail/alps.hpp"

#ifdef USE_MTM
#include "types/mt_matrix/mt_matrix.h"
#include "types/mt_matrix/algorithms.hpp"
#endif

#include "dmrg/models/factory.h"
#include "dmrg/models/coded/factory.h"
#include "dmrg/models/continuum/factory.h"
#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#include "dmrg/models/alps/model.hpp"
#endif

#include <fstream>
#include <boost/tokenizer.hpp>

// BLAS matrix
typedef alps::numeric::matrix<double> matrix;
typedef alps::numeric::matrix<std::complex<double> > cmatrix;

// MT matrix
#ifdef USE_MTM
typedef alps::numeric::mt_matrix<double> mtmatrix;
typedef alps::numeric::mt_matrix<std::complex<double> > cmtmatrix;
#endif

// parallel matrix
#ifdef USE_AMBIENT
typedef ambient::numeric::matrix<double> pmatrix;
typedef ambient::numeric::matrix<std::complex<double> > cpmatrix;
#endif

// Definition of init function
template <class Matrix, class SymmGroup>
void init_model_parser();

// init MACROS
#define impl_init_model(MATRIX, SYMMGROUP)														\
template <>																						\
void init_model_parser<MATRIX, SYMMGROUP> ()													\
{																								\
    Lattice_ptr lat;																			\
    model_traits<MATRIX, SYMMGROUP>::model_ptr phys_model;                                      \
	SYMMGROUP::charge initc;																	\
	BaseParameters parms;																		\
	model_parser<MATRIX,SYMMGROUP>("", "", parms, lat, phys_model);                             \
}



// Implementations

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
template <class Matrix, class SymmGroup>
void model_parser (std::string lattice_lib, std::string model_lib,
                   BaseParameters & parms,
                   Lattice_ptr & lattice,
                   typename model_traits<Matrix, SymmGroup>::model_ptr & model)
{
    // Lattice
    if (lattice_lib == "alps") {
#ifdef ENABLE_ALPS_MODELS
        lattice = Lattice_ptr(new ALPSLattice(parms));
#else
        throw std::runtime_error("This code was compiled without alps lattice.");
#endif
    } else if (lattice_lib == "coded") {
        lattice = lattice_factory(parms);
    } else if (lattice_lib == "continuum") {
        lattice = cont_lattice_factory(parms);
    } else {
        throw std::runtime_error("Don't know this lattice_library!");
    }

    // Model
    if (model_lib == "alps") {
#ifdef ENABLE_ALPS_MODELS
        if (lattice_lib != "alps")
            throw std::runtime_error("ALPS models require ALPS lattice.");
        model = typename model_traits<Matrix, SymmGroup>::model_ptr(
                    new ALPSModel<Matrix, SymmGroup>(static_cast<ALPSLattice*>(lattice.get())->alps_graph(),
                                                     parms)
                );
#else
        throw std::runtime_error("This code was compiled without alps models.");
#endif
    } else if (model_lib == "coded") {
        model = model_factory<Matrix, SymmGroup>::parse(*lattice, parms);
    } else if (model_lib == "continuum") {
        model = cont_model_factory<Matrix, SymmGroup>::parse(*lattice, parms);
    } else {
        throw std::runtime_error("Don't know this model_library!");
    }
    
}

