/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "block_matrix/block_matrix.h"
#include "block_matrix/block_matrix_algorithms.h"
#include "block_matrix/symmetry.h"

#include <vector>
#include <string>
#include <iostream>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#include <alps/utility/encode.hpp>
#endif

#include "lattice.h"

#include "generate_mpo.hpp"

namespace app {
    
    template<class Matrix, class SymmGroup>
    struct Measurement_Term
    {
        typedef block_matrix<Matrix, SymmGroup> op_t;
        enum type_t {Local, Average, Correlation, HalfCorrelation, CorrelationNN, HalfCorrelationNN};
        
        std::string name;
        type_t type;
        std::vector<std::pair<op_t, bool> > operators;
        op_t fill_operator;
    };
    
    // implement me for your purpose!
    template<class Matrix, class SymmGroup>
    class Measurements
    {
    public:
        int n_terms() const
        {
        	return terms.size();
        }
        const Measurement_Term<Matrix, SymmGroup>& operator[](int i) const
        {
        	return terms[i];
        }
        
        const typename Measurement_Term<Matrix, SymmGroup>::op_t& get_identity() const
        {
        	return ident;
        }
        
        void set_identity(const typename Measurement_Term<Matrix, SymmGroup>::op_t& ident_)
        {
        	ident = ident_;
        }
        
        void add_term (Measurement_Term<Matrix, SymmGroup> const & term)
        {
        	terms.push_back(term);
        }
        
    protected:
        std::vector<Measurement_Term<Matrix, SymmGroup> > terms;
        typename Measurement_Term<Matrix, SymmGroup>::op_t ident;
    };
}

#include "meas_detail.hpp"

// call me to do something!

namespace app {
    
    
    template<class Matrix, class SymmGroup>
    void measure(MPS<Matrix, SymmGroup> & mps, Lattice const & lat,
    			 Measurements<Matrix, SymmGroup> const & meas, alps::hdf5::oarchive & ar)
    {
    	std::string basepath = "/spectrum/results/";
    	for (int i = 0; i < meas.n_terms(); ++i)
    	{
    		switch (meas[i].type)
    		{
                case Measurement_Term<Matrix, SymmGroup>::Local:
                    assert(meas[i].operators.size() == 1);
                    meas_detail::measure_local(mps, lat,
                                               meas.get_identity(), meas[i].operators[0].first,
                                               ar, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Average:
                    assert(meas[i].operators.size() == 1);
                    meas_detail::measure_average(mps, lat,
                                                 meas.get_identity(), meas[i].operators[0].first,
                                                 ar, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Correlation:
                    meas_detail::measure_correlation<generate_mpo::CorrMaker<Matrix, SymmGroup> >(mps, lat, meas.get_identity(),
                                                                                                  meas[i].fill_operator, meas[i].operators,
                                                                                                  ar, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
                    meas_detail::measure_correlation<generate_mpo::CorrMaker<Matrix, SymmGroup> >(mps, lat, meas.get_identity(),
                                                                                                  meas[i].fill_operator, meas[i].operators,
                                                                                                  ar, basepath + alps::hdf5_name_encode(meas[i].name), true);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
#ifndef NDEBUG
                    if (meas[i].operators.size() % 2 != 0)
                        std::runtime_error("Next neighbors correlators have to have even number of operators");
#endif
                    meas_detail::measure_correlation<generate_mpo::CorrMakerNN<Matrix, SymmGroup> >(mps, lat, meas.get_identity(),
                                                                                                    meas[i].fill_operator, meas[i].operators,
                                                                                                    ar, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
#ifndef NDEBUG
                    if (meas[i].operators.size() % 2 != 0)
                        std::runtime_error("Next neighbors correlators have to have even number of operators");
#endif
                    meas_detail::measure_correlation<generate_mpo::CorrMakerNN<Matrix, SymmGroup> >(mps, lat, meas.get_identity(),
                                                                                                    meas[i].fill_operator, meas[i].operators,
                                                                                                    ar, basepath + alps::hdf5_name_encode(meas[i].name), true);
                    break;
    		}
    	}
    }
    
} // namespace


// Outout operators

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, typename app::Measurement_Term<Matrix, SymmGroup>::type_t const & t)
{
	switch (t)
	{
        case app::Measurement_Term<Matrix, SymmGroup>::Local:
            os << "Local";
            break;
        case app::Measurement_Term<Matrix, SymmGroup>::Average:
            os << "Average";
            break;
        case app::Measurement_Term<Matrix, SymmGroup>::Correlation:
            os << "Correlation";
            break;
        case app::Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
            os << "HalfCorrelation";
            break;
        case app::Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
            os << "CorrelationNN";
            break;
        case app::Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
            os << "HalfCorrelationNN";
            break;
	}
	return os;
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, std::vector<std::pair<block_matrix<Matrix, SymmGroup>, bool> > const & v)
{
	for (int i=0; i<v.size(); ++i) {
		std::string t = (v[i].second) ? "fermionic " : "";
		os << t << "operator" << std::endl;
		os << v[i].first;
	}
	return os;
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, app::Measurement_Term<Matrix, SymmGroup> const & term)
{
	os << "** MEASUREMENT: " << term.name << " **" << std::endl;
	os << " - type: " << term.type << std::endl;
	os << " - Fill operator:" << std::endl << term.operators;
	return os;
}
template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, app::Measurements<Matrix, SymmGroup> const & meas)
{
	os << "** IDENTIY **" << std::endl;
	os << meas.get_identity();
	for (int i=0; i<meas.n_terms(); ++i) {
		os << meas[i] << std::endl;
	}
	return os;
}



#endif
