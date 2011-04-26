
#define NO_ZOUT_IN_HEADERS

#include <complex>
#include <vector>
#include "dense_matrix/dense_matrix.h"
#include "dense_matrix/matrix_interface.hpp"
#include "dense_matrix/resizable_matrix_interface.hpp"
#include "dense_matrix/dense_matrix_algorithms.h"
#include "dense_matrix/matrix_algorithms.hpp"
#include "dense_matrix/dense_matrix_blas.hpp"
#include "dense_matrix/aligned_allocator.h"

#include "models.h"

#include <fstream>
#include <boost/tokenizer.hpp>

#include "hamiltonians.hpp"
#include "alps_lattice.hpp"
#include "alps_model.hpp"

#include "utils/DmrgParameters.h"
#include <alps/parameter.h>



//Lattice * lattice_factory (std::string const & lattice, std::ifstream & ifs)
//{
//    if (lattice == std::string("alps_lattice"))
//        return new ALPSLattice(ifs);
//    else {
//        throw std::runtime_error("Don't know this lattice!");
//        return NULL;
//    }
//}

namespace app {

#define model_parser_for(MATRIX, SYMMGROUP) 							\
	{																	\
		Hamiltonian<MATRIX, SYMMGROUP>* H;								\
		SYMMGROUP::charge initc;										\
		model_parser(std::string(), std::string(), lat, H, initc);		\
	}
void init_model_parser()
{
	Lattice * lat;

	typedef blas::dense_matrix<double> matrix1;
	typedef blas::dense_matrix<double, std::vector<double, aligned_allocator<double> > > matrix2;
	typedef blas::dense_matrix<std::complex<double> > cmatrix1;
	typedef blas::dense_matrix<std::complex<double>, std::vector<std::complex<double>, aligned_allocator<std::complex<double> > > > cmatrix2;

	model_parser_for(matrix1, U1)
	model_parser_for(matrix1, TwoU1)
	model_parser_for(matrix2, U1)
	model_parser_for(matrix2, TwoU1)
	model_parser_for(cmatrix1, U1)
	model_parser_for(cmatrix1, TwoU1)
	model_parser_for(cmatrix2, U1)
	model_parser_for(cmatrix2, TwoU1)

}
#undef model_parser_for

// Implementations

typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

template<class Matrix, class SymmGroup>
struct hamil_factory {
	static Hamiltonian<Matrix, SymmGroup> * parse (BaseParameters & model, Lattice const & lattice);
};

template <class SymmGroup>
typename SymmGroup::charge init_qn (ModelParameters & model);


template <class Matrix, class SymmGroup>
void model_parser (std::string const & type, std::string const & fname,
                   Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H, typename SymmGroup::charge& initc)
{

    if (type == "alps") {
        std::ifstream ifs(fname.c_str());
        alps::Parameters parms(ifs);
        ALPSLattice * lat_tmp = new ALPSLattice(parms);
        ALPSModel<Matrix, SymmGroup> * Htmp = new ALPSModel<Matrix, SymmGroup>(lat_tmp->alps_graph(), parms);
        lattice = lat_tmp;
        H = Htmp;
        initc = Htmp->init_qn(parms);
    } else if (type == "coded") {
        std::ifstream model_file(fname.c_str());
        std::ifstream ifs(fname.c_str());
        alps::Parameters parms(ifs);
        lattice = new ALPSLattice(parms);
        ModelParameters model(model_file);
        H = hamil_factory<Matrix, SymmGroup>::parse(model, *lattice);
        initc = init_qn<SymmGroup>(model);
    } else {
        throw std::runtime_error("Don't know this type of lattice / model!");
    }
    
}


template<class Matrix>
struct hamil_factory<Matrix, TwoU1> {
	static Hamiltonian<Matrix, TwoU1> * parse (BaseParameters & model, Lattice const & lattice)
	{
		if (model.get<std::string>("model") == std::string("fermi_hubbard"))
			return new TwoU1_FermiHubbard<Matrix>(lattice, model);
		else {
			throw std::runtime_error("Don't know this model!");
			return NULL;
		}
	}
};

template<class Matrix>
struct hamil_factory<Matrix, U1> {
	static Hamiltonian<Matrix, U1> * parse (BaseParameters & model, Lattice const & lattice)
	{
		if (model.get<std::string>("model") == std::string("heisenberg"))
			return new Heisenberg<Matrix>(lattice, model.get<double>("Jxy"), model.get<double>("Jz"));
		else if (model.get<std::string>("model") == std::string("HCB"))
			return new HCB<Matrix>(lattice);
		else if (model.get<std::string>("model") == std::string("FreeFermions"))
			return new FreeFermions<Matrix>(lattice, model.get<double>("t"));
		else {
			throw std::runtime_error("Don't know this model!");
			return NULL;
		}
	}
};


template <>
U1::charge init_qn<U1> (ModelParameters & model)
{
    return model.get<int>("u1_total_charge");
}
template <>
TwoU1::charge init_qn<TwoU1> (ModelParameters & model)
{
	TwoU1::charge initc;
    initc[0] = model.get<int>("u1_total_charge1");
    initc[1] = model.get<int>("u1_total_charge2");
	return initc;
}

    
} // namespace
