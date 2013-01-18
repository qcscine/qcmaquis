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

#include <boost/regex.hpp>
#include <boost/static_assert.hpp>


template<class Matrix, class SymmGroup>
class Measurement_Term
{
public:
    typedef block_matrix<Matrix, SymmGroup> op_t;
    enum type_t {Local, MPSBonds, Average, Correlation, HalfCorrelation, CorrelationNN, HalfCorrelationNN, Custom, Overlap};

    Measurement_Term() {}

    std::string name;
    type_t type;
    std::vector<std::pair<op_t, bool> > operators;
    op_t fill_operator;
    
    // this is somewhat unlucky interface-wise, to say the least
    std::vector< std::vector< std::pair<int, op_t> > > custom_ops;
    
    Measurement_Term<Matrix, SymmGroup> * clone() const
    {
        return do_clone();
    }
protected:
    virtual Measurement_Term<Matrix, SymmGroup> * do_clone() const
    {
        return new Measurement_Term<Matrix, SymmGroup>(*this);
    }
    
    Measurement_Term(Measurement_Term const & m)
    : name(m.name)
    , type(m.type)
    , operators(m.operators)
    , fill_operator(m.fill_operator)
    { }
};
template<class Matrix, class SymmGroup>
inline Measurement_Term<Matrix, SymmGroup>* new_clone( const Measurement_Term<Matrix, SymmGroup>& m )
{
    return m.clone();
}


template <class Matrix, class SymmGroup>
class OverlapMeasurement : public Measurement_Term<Matrix, SymmGroup> {
public:
    //    typedef typename Measurement_Term<Matrix, SymmGroup>::op_t op_t;
    typedef Measurement_Term<Matrix, SymmGroup> base;
    
    OverlapMeasurement()
    { this->type = Measurement_Term<Matrix, SymmGroup>::Overlap; }
    
    // checkpoint of the <bra| wavefunction
    std::string bra_ckp;
protected:
    virtual Measurement_Term<Matrix, SymmGroup> * do_clone() const
    {
        return new OverlapMeasurement<Matrix, SymmGroup>(*this);
    }
};


// implement me for your purpose!
template<class Matrix, class SymmGroup>
class Measurements
{
public:
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    typedef typename mterm_t::type_t meas_type_t;
    typedef typename mterm_t::op_t op_t;
    
    enum system_type_t {Wavefunction, Densitymatrix};
    
    Measurements(system_type_t system_type=Wavefunction) : super_meas( system_type == Densitymatrix ) { }
    
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
    	terms.push_back(term.clone());
    }

    void clear ()
    {
    	terms.clear();
        ident = op_t();
    }
    
    bool is_super_meas() const
    {
        return super_meas;
    }
    
protected:
    boost::ptr_vector<mterm_t> terms;
    op_t ident;
    bool super_meas;
};


#include "meas_detail.hpp"

template<class Matrix, class SymmGroup>
void measure_on_mps(MPS<Matrix, SymmGroup> & mps, Lattice const & lat,
                    Measurements<Matrix, SymmGroup> const & meas,
                    std::string const & h5name, std::string basepath = std::string("/spectrum/results/"))
{
	
    if (meas.n_terms() > 0) {
        bool super_meas=meas.is_super_meas();
        
        boost::scoped_ptr<meas_detail::LocalMPSMeasurement<Matrix, SymmGroup> > local_measurement;
        if (!super_meas)
            local_measurement.reset( new meas_detail::LocalMPSMeasurement<Matrix, SymmGroup>(mps, lat) );
        
        for (int i = 0; i < meas.n_terms(); ++i)
        {
            maquis::cout << "Calculating " << meas[i].name << std::endl;
            switch (meas[i].type)
            {
                case Measurement_Term<Matrix, SymmGroup>::Local:
                    assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                    if (!super_meas && meas[i].operators.size() == 1) // Local measurements are fast and efficient!
                        local_measurement->site_term(meas[i].operators[0],
                                                    h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                    else
                        meas_detail::measure_local(mps, lat,
                                                   meas.get_identity(), meas[i].fill_operator,
                                                   meas[i].operators,
                                                   h5name, basepath + alps::hdf5_name_encode(meas[i].name), super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::MPSBonds:
                    assert(meas[i].operators.size() == 2);
                    local_measurement->bond_term(meas[i].operators,
                                                h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Average:
                    assert(meas[i].operators.size() == 1  || meas[i].operators.size() == 2);
                    meas_detail::measure_average(mps, lat,
                                                 meas.get_identity(), meas[i].fill_operator,
                                                 meas[i].operators,
                                                 h5name, basepath + alps::hdf5_name_encode(meas[i].name), super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Overlap:
                    meas_detail::measure_overlap(mps, dynamic_cast<OverlapMeasurement<Matrix, SymmGroup> const & >(meas[i]).bra_ckp,
                                                 h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Correlation:
                    meas_detail::measure_correlation(mps, lat, meas.get_identity(),
                                                     meas[i].fill_operator, meas[i].operators,
                                                     h5name, basepath + alps::hdf5_name_encode(meas[i].name), false, false, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelation:
                    meas_detail::measure_correlation(mps, lat, meas.get_identity(),
                                                     meas[i].fill_operator, meas[i].operators,
                                                     h5name, basepath + alps::hdf5_name_encode(meas[i].name), true, false, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::CorrelationNN:
                    if (meas[i].operators.size() % 2 != 0)
                        throw std::runtime_error("Next neighbors correlators have to have even number of operators");
                    meas_detail::measure_correlation(mps, lat, meas.get_identity(),
                                                     meas[i].fill_operator, meas[i].operators,
                                                     h5name, basepath + alps::hdf5_name_encode(meas[i].name), false, true, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::HalfCorrelationNN:
                    if (meas[i].operators.size() % 2 != 0)
                        throw std::runtime_error("Next neighbors correlators have to have even number of operators");
                    meas_detail::measure_correlation(mps, lat, meas.get_identity(),
                                                     meas[i].fill_operator, meas[i].operators,
                                                     h5name, basepath + alps::hdf5_name_encode(meas[i].name), true, true, super_meas);
                    break;
                case Measurement_Term<Matrix, SymmGroup>::Custom:
                    meas_detail::measure_custom(mps, lat, meas.get_identity(),
                                                meas[i].fill_operator, meas[i].custom_ops,
                                                h5name, basepath + alps::hdf5_name_encode(meas[i].name));
                    break;
            }
        }
    }
}

template <class Matrix, class SymmGroup>
void parse_overlaps(BaseParameters const & parms, size_t sweep, Measurements<Matrix, SymmGroup> & meas)
{
    /* Syntax for MEASURE_OVERLAP:
     *  (1) MEASURE_OVERLAP[obsname] = "/path/to/ckp.h5"
     *  (2) MEASURE_OVERLAP[obsname(sweep)] = "/path/to/ckp.h5"
     *
     * `obsname` is the name that will be given in hdf5 output.
     * if `sweep` is prensent, the overlap will only be computed when the sweep number
     * matches the given one.
     */
    boost::regex expression("^MEASURE_OVERLAP\\[([a-zA-Z]+)(\\(([0-9]+)\\))?\\]$");
    boost::smatch what;
    for (BaseParameters::const_iterator it=parms.begin();it != parms.end();++it) {
        std::string lhs = it->key();
        if (boost::regex_match(lhs, what, expression)) {
            if (!what[3].matched || boost::lexical_cast<long>(what.str(3)) != sweep)
                continue;
            
            OverlapMeasurement<Matrix, SymmGroup> term;
            term.name = what.str(1);
            term.bra_ckp = it->value();
            meas.add_term(term);
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
        case Measurement_Term<Matrix, SymmGroup>::Overlap:
            os << "Overlap";
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
std::ostream& operator<<(std::ostream& os, OverlapMeasurement<Matrix, SymmGroup> const & term)
{
	os << "** MEASUREMENT: " << term.name << " **" << std::endl;
	os << " - type: " << term.type << std::endl;
	os << " - Path to |bra>:" << term.bra_ckp << std::endl;
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
