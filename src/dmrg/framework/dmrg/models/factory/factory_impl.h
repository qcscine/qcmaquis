/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include <complex>
#include <vector>
#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/resizable_matrix_interface.hpp"
#include "types/dense_matrix/dense_matrix_algorithms.h"
#include "types/dense_matrix/matrix_algorithms.hpp"
#include "types/dense_matrix/algorithms.hpp"
#include "types/dense_matrix/dense_matrix_blas.hpp"
#include "types/dense_matrix/aligned_allocator.h"

#ifdef USE_MTM
#include "types/dense_matrix/mt_matrix.h"
#endif

#include "dmrg/models/factory.h"
#include "dmrg/models/coded/lattice.hpp"
#include "dmrg/models/coded/hamiltonians.hpp"
#include "dmrg/models/coded/measurements.hpp"
#ifdef ENABLE_ALPS_MODELS
#include "dmrg/models/alps/lattice.hpp"
#include "dmrg/models/alps/model.hpp"
#endif

#include <fstream>
#include <boost/tokenizer.hpp>

namespace app {

	// BLAS Matrix
	typedef maquis::types::dense_matrix<double> matrix1;
	typedef maquis::types::dense_matrix<std::complex<double> > cmatrix1;
	
	// MT Matrix
#ifdef USE_MTM
	typedef mt_matrix<double> mtmatrix1;
	typedef mt_matrix<std::complex<double> > cmtmatrix1;
#endif

	// Definition of init function
	template <class Matrix, class SymmGroup>
    void init_model_parser();

// init MACROS
#define impl_init_model(MATRIX, SYMMGROUP)															\
	template <>																						\
	void init_model_parser<MATRIX, SYMMGROUP> ()													\
	{																								\
        Lattice * lat;																				\
		Hamiltonian<MATRIX, SYMMGROUP> H;															\
		Measurements<MATRIX, SYMMGROUP> meas;														\
		SYMMGROUP::charge initc;																	\
		BaseParameters parms;																		\
		model_parser("", "", parms, lat, H, initc, meas);                                           \
	}
    
    
    
    // Implementations
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    template<class Matrix, class SymmGroup>
    struct hamil_factory {
        static Hamiltonian<Matrix, SymmGroup> parse (BaseParameters & model, Lattice const & lattice);
    };
    inline Lattice * lattice_factory (BaseParameters & model);
    
    template <class SymmGroup>
    typename SymmGroup::charge init_qn (BaseParameters & model);
    
    template <class Matrix, class SymmGroup>
    void model_parser (std::string lattice_lib, std::string model_lib,
                       BaseParameters & parms,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>& H,
                       typename SymmGroup::charge& initc, Measurements<Matrix, SymmGroup>& meas)
    {
        // Lattice
        if (lattice_lib == "alps") {
#ifdef ENABLE_ALPS_MODELS
            lattice = new ALPSLattice(parms);
#else
            throw std::runtime_error("This code was compiled without alps lattice.");
#endif
        } else if (lattice_lib == "coded") {
            lattice = lattice_factory(parms);
        } else {
            throw std::runtime_error("Don't know this lattice_library!");
        }

        // Model
        if (model_lib == "alps") {
#ifdef ENABLE_ALPS_MODELS
            if (lattice_lib != "alps")
                throw std::runtime_error("ALPS models require ALPS lattice.");
            
            ALPSModel<Matrix, SymmGroup> HALPS = ALPSModel<Matrix, SymmGroup>(static_cast<ALPSLattice*>(lattice)->alps_graph(),
                                                                              parms);
            H = HALPS.get_hamiltonian();
            initc = HALPS.init_qn(parms);
            meas = HALPS.parse_measurements(parms);
#else
            throw std::runtime_error("This code was compiled without alps models.");
#endif
        } else if (model_lib == "coded") {
            H = hamil_factory<Matrix, SymmGroup>::parse(parms, *lattice);
            initc = init_qn<SymmGroup>(parms);
            meas = CodedMeasurements<Matrix, SymmGroup>(*lattice, parms);
        } else {
            throw std::runtime_error("Don't know this model_library!");
        }
        
    }
    
    
    
    inline Lattice * lattice_factory (BaseParameters & model)
    {
        if (model.get<std::string>("LATTICE") == std::string("periodic chain lattice"))
            return new ChainLattice(model, true);
        else if (model.get<std::string>("LATTICE") == std::string("open chain lattice"))
            return new ChainLattice(model, false);
        else {
            throw std::runtime_error("Don't know this lattice!");
        }
    }

    
} // namespace
