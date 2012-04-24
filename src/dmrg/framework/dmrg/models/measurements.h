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

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include <vector>
#include <string>
#include <iostream>

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#include <alps/utility/encode.hpp>
#endif

#include "dmrg/models/lattice.h"
#include "dmrg/models/generate_mpo.hpp"

#include <stdexcept>

template<class Matrix, class SymmGroup>
struct Measurement_Term
{
    typedef block_matrix<Matrix, SymmGroup> op_t;
    enum type_t {Local, MPSBonds, Average, Correlation, HalfCorrelation, CorrelationNN, HalfCorrelationNN};
    
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
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    typedef typename mterm_t::type_t meas_type_t;
    typedef typename mterm_t::op_t op_t;
            
    int n_terms() const
    {
    	return terms.size();
    }
    const mterm_t& operator[](int i) const
    {
    	return terms[i];
    }
    const mterm_t& get(std::string const & name) const
    {
        for (int i=0; i<terms.size(); ++i) {
            if (terms[i].name == name)
                return terms[i];
        }
        throw std::runtime_error("Measurement "+name+" not found!");
        return *(terms.end());
    }
    
    const op_t& get_identity() const
    {
    	return ident;
    }
    
    void set_identity(const op_t& ident_)
    {
    	ident = ident_;
    }
    
    void add_term (mterm_t const & term)
    {
    	terms.push_back(term);
    }

    void clear ()
    {
    	terms.clear();
        ident = op_t();
    }

protected:
    std::vector<mterm_t> terms;
    op_t ident;
};


#include "meas_detail.hpp"

template<class Matrix, class SymmGroup>
void measure_on_mps(MPS<Matrix, SymmGroup> & mps, Lattice const & lat,
			 Measurements<Matrix, SymmGroup> const & meas,
             std::string const & h5name, std::string basepath = std::string("/spectrum/results/"))
{
	
    meas_detail::LocalMPSMeasurement<Matrix, SymmGroup> local_measurement(mps, lat);
    
    for (int i = 0; i < meas.n_terms(); ++i)
	{
		maquis::cout << "Calculating " << meas[i].name << std::endl;
		switch (meas[i].type)
		{
            case Measurement_Term<Matrix, SymmGroup>::Local:
                assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                if (meas[i].operators.size() == 1) // Local measurements are fast and efficient!
                    local_measurement.site_term(meas[i].operators[0],
                                                h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                else
                    meas_detail::measure_local(mps, lat,
                                               meas.get_identity(), meas[i].fill_operator,
                                               meas[i].operators,
                                               h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                break;
            case Measurement_Term<Matrix, SymmGroup>::MPSBonds:
                assert(meas[i].operators.size() == 2);
                local_measurement.bond_term(meas[i].operators,
                                            h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                break;
            case Measurement_Term<Matrix, SymmGroup>::Average:
                assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                meas_detail::measure_average(mps, lat,
                                             meas.get_identity(), meas[i].fill_operator,
                                             meas[i].operators,
                                             h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                break;
            case Measurement_Term<Matrix, SymmGroup>::Correlation:
                meas_detail::measure_correlation(mps, lat, meas.get_identity(),
						  	  	  	  	  	  	 meas[i].fill_operator, meas[i].operators,
											  	 h5name, basepath + alps::hdf5_name_encode(meas[i].name), false, false);
                break;
            case Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
                meas_detail::measure_correlation(mps, lat, meas.get_identity(),
												  meas[i].fill_operator, meas[i].operators,
												  h5name, basepath + alps::hdf5_name_encode(meas[i].name), true, false);
                break;
            case Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
#ifndef NDEBUG
                if (meas[i].operators.size() % 2 != 0)
                    throw std::runtime_error("Next neighbors correlators have to have even number of operators");
#endif
                meas_detail::measure_correlation(mps, lat, meas.get_identity(),
												meas[i].fill_operator, meas[i].operators,
												h5name, basepath + alps::hdf5_name_encode(meas[i].name), false, true);
                break;
            case Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
#ifndef NDEBUG
                if (meas[i].operators.size() % 2 != 0)
                    throw std::runtime_error("Next neighbors correlators have to have even number of operators");
#endif
                meas_detail::measure_correlation(mps, lat, meas.get_identity(),
												meas[i].fill_operator, meas[i].operators,
												h5name, basepath + alps::hdf5_name_encode(meas[i].name), true, true);
                break;
		}
	}
}

// Output operators

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, typename Measurement_Term<Matrix, SymmGroup>::type_t const & t)
{
	switch (t)
	{
        case Measurement_Term<Matrix, SymmGroup>::Local:
            os << "Local";
            break;
        case Measurement_Term<Matrix, SymmGroup>::Average:
            os << "Average";
            break;
        case Measurement_Term<Matrix, SymmGroup>::Correlation:
            os << "Correlation";
            break;
        case Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
            os << "HalfCorrelation";
            break;
        case Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
            os << "CorrelationNN";
            break;
        case Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
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
std::ostream& operator<<(std::ostream& os, Measurement_Term<Matrix, SymmGroup> const & term)
{
	os << "** MEASUREMENT: " << term.name << " **" << std::endl;
	os << " - type: " << term.type << std::endl;
	os << " - Fill operator:" << std::endl << term.operators;
	return os;
}

template<class Matrix, class SymmGroup>
std::ostream& operator<<(std::ostream& os, Measurements<Matrix, SymmGroup> const & meas)
{
	os << "** IDENTIY **" << std::endl;
	os << meas.get_identity();
	for (int i=0; i<meas.n_terms(); ++i) {
		os << meas[i] << std::endl;
	}
	return os;
}

#endif
