/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#include <complex>
#include <vector>
#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/aligned_allocator.h"

#ifdef USE_MTM
#include "dense_matrix/mt_matrix.h"
#endif

#include "app/models.h"
#include "app/hamiltonians.hpp"
#include "app/measurements.hpp"
#include "app/alps_lattice.hpp"
#include "app/alps_model.hpp"

#include <alps/parameter.h>
#include <alps/expression.h>

#include <fstream>
#include <boost/tokenizer.hpp>

namespace app {

	// BLAS Matrix
	typedef blas::dense_matrix<double> matrix1;
	typedef blas::dense_matrix<std::complex<double> > cmatrix1;

	// BLAS Matrix with aligned allocator
	typedef blas::dense_matrix<double, std::vector<double, aligned_allocator<double> > > matrix2;
	typedef blas::dense_matrix<std::complex<double>, std::vector<std::complex<double>, aligned_allocator<std::complex<double> > > > cmatrix2;
	
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
		model_parser(std::string(), parms, lat, H, initc, meas);									\
	}
    
    
    
    // Implementations
    
    typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
    
    template<class Matrix, class SymmGroup>
    struct hamil_factory {
        static Hamiltonian<Matrix, SymmGroup> parse (BaseParameters & model, Lattice const & lattice);
    };
    
    template <class SymmGroup>
    typename SymmGroup::charge init_qn (BaseParameters & model);
    
    template <class Matrix, class SymmGroup>
    void model_parser (std::string const & type, BaseParameters & parms,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>& H,
                       typename SymmGroup::charge& initc, Measurements<Matrix, SymmGroup>& meas)
    {
        
        if (type == "alps") {
            ALPSLattice * lat_tmp = new ALPSLattice(parms);
            ALPSModel<Matrix, SymmGroup> HALPS = ALPSModel<Matrix, SymmGroup>(lat_tmp->alps_graph(), parms);
            lattice = lat_tmp;
            H = HALPS.get_hamiltonian();
            initc = HALPS.init_qn(parms);
            meas = HALPS.parse_measurements(parms);
        } else if (type == "coded") {
            lattice = new ALPSLattice(parms);
            H = hamil_factory<Matrix, SymmGroup>::parse(parms, *lattice);
            initc = init_qn<SymmGroup>(parms);
            meas = CodedMeasurements<Matrix, SymmGroup>(*lattice, parms);
        } else {
            throw std::runtime_error("Don't know this type of lattice / model!");
        }
        
    }
    
    
} // namespace
