
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

#include "dense_matrix/mt_matrix.h"

#include "models.h"

#include <fstream>
#include <boost/tokenizer.hpp>

#include "hamiltonians.hpp"
#include "measurements.hpp"
#include "alps_lattice.hpp"
#include "alps_model.hpp"

#include "utils/DmrgParameters.h"
#include <alps/parameter.h>
#include <alps/expression.h>


namespace app {

#define model_parser_for(MATRIX, SYMMGROUP) 														\
	{																								\
		Hamiltonian<MATRIX, SYMMGROUP>* H;															\
		Measurements<MATRIX, SYMMGROUP> meas;														\
		SYMMGROUP::charge initc;																	\
		ModelParameters model = model_parser(std::string(), std::string(), lat, H, initc, meas);	\
	}
    void init_model_parser()
    {
        Lattice * lat;
        
        // BLAS Matrix
        typedef blas::dense_matrix<double> matrix1;
        typedef blas::dense_matrix<std::complex<double> > cmatrix1;

        // BLAS Matrix with aligned allocator
        typedef blas::dense_matrix<double, std::vector<double, aligned_allocator<double> > > matrix2;
        typedef blas::dense_matrix<std::complex<double>, std::vector<std::complex<double>, aligned_allocator<std::complex<double> > > > cmatrix2;
        
        // MT Matrix
        typedef mt_matrix<double> mtmatrix1;
        typedef mt_matrix<std::complex<double> > cmtmatrix1;
        
        
#ifdef USE_MTM
        model_parser_for(mtmatrix1, U1)
        model_parser_for(mtmatrix1, TwoU1)        
        //        model_parser_for(cmtmatrix1, U1)
        //        model_parser_for(cmtmatrix1, TwoU1)
#else
        model_parser_for(matrix1, U1)
        model_parser_for(matrix1, TwoU1)
        model_parser_for(matrix2, U1)
        model_parser_for(matrix2, TwoU1)
        model_parser_for(cmatrix1, U1)
        model_parser_for(cmatrix1, TwoU1)
        model_parser_for(cmatrix2, U1)
        model_parser_for(cmatrix2, TwoU1)
#endif
        
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
    
    ModelParameters convert_parms (alps::Parameters const & parms);
    
    template <class Matrix, class SymmGroup>
    ModelParameters model_parser (std::string const & type, std::string const & fname,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H,
                       typename SymmGroup::charge& initc, Measurements<Matrix, SymmGroup>& meas)
    {
        
        if (type == "alps") {
            std::ifstream ifs(fname.c_str());
            alps::Parameters parms(ifs);
            ALPSLattice * lat_tmp = new ALPSLattice(parms);
            ALPSModel<Matrix, SymmGroup> * Htmp = new ALPSModel<Matrix, SymmGroup>(lat_tmp->alps_graph(), parms);
            lattice = lat_tmp;
            H = Htmp;
            initc = Htmp->init_qn(parms);
            meas = Htmp->parse_measurements(parms);
            return convert_parms(parms);
        } else if (type == "coded") {
            std::ifstream model_file(fname.c_str());
            std::ifstream ifs(fname.c_str());
            alps::Parameters parms(ifs);
            lattice = new ALPSLattice(parms);
            ModelParameters model = convert_parms(parms);
            H = hamil_factory<Matrix, SymmGroup>::parse(model, *lattice);
            initc = init_qn<SymmGroup>(model);
            meas = CodedMeasurements<Matrix, SymmGroup>(*lattice, model);
            return model;
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
    
    
    ModelParameters convert_parms (alps::Parameters const & parms)
    {
    	ModelParameters model;
        alps::expression::ParameterEvaluator<double> eval(parms,false);
        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            try {
                alps::expression::Expression<double> expr(it->value());
                if (expr.can_evaluate(eval)) {
                    double value = expr.value(eval);
                    if (alps::numeric::is_zero(value - static_cast<double>(static_cast<int>(value))))
                    	model.set(it->key(), static_cast<int>(value+ (value > 0 ? 0.25 : -0.25)));
                    else
                    	model.set(it->key(), value);
                } else {
                    expr.partial_evaluate(eval);
                    model.set(it->key(), boost::lexical_cast<std::string>(expr));
                }
            } catch(...) {
              // we had a problem evaluating, use original full value
            	model.set(it->key(), boost::lexical_cast<std::string>(it->value()));
            }
        }
        return model;
    }


} // namespace
