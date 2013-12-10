/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MEASUREMENTS_H
#define MEASUREMENTS_H

#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"

#include <vector>
#include <string>
#include <iostream>

#include <stdexcept>

#include <boost/regex.hpp>
#include <boost/static_assert.hpp>


template<class Matrix, class SymmGroup>
class Measurement_Term
{
public:
    typedef block_matrix<Matrix, SymmGroup> op_t;
    enum type_t {Local, LocalAt, MPSBonds, Average, Correlation, HalfCorrelation, CorrelationNN, HalfCorrelationNN, Custom, Overlap, DMOverlap, DMMultioverlap};

    Measurement_Term() {}

    std::string name;
    type_t type;
    std::vector<std::pair<std::vector<op_t>, bool> > operators;
    
    // this is somewhat unlucky interface-wise, to say the least
    // Custom: all inner vector are summed together
    std::vector< std::vector< std::pair<int, op_t> > > custom_ops;
    
    // used by LocalAt for the positions where ops are evaluated
    std::vector< std::vector<std::size_t> > positions;
    
    // physical_i for wavefunction in case of density matrix measurements
    Index<SymmGroup> phys_psi;
    
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
    , custom_ops(m.custom_ops)
    , positions(m.positions)
    , phys_psi(m.phys_psi)
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

template <class Matrix, class SymmGroup>
class DMOverlapMeasurement : public Measurement_Term<Matrix, SymmGroup> {
public:
    //    typedef typename Measurement_Term<Matrix, SymmGroup>::op_t op_t;
    typedef Measurement_Term<Matrix, SymmGroup> base;
    
    MPS<Matrix, SymmGroup> mps_ident;
    std::vector<MPS<Matrix, SymmGroup> > overlaps_mps;
    std::vector<std::string> labels;
    
protected:
    virtual Measurement_Term<Matrix, SymmGroup> * do_clone() const
    {
        return new DMOverlapMeasurement<Matrix, SymmGroup>(*this);
    }
};

namespace detail {
    class name_not_in_list {
    public:
        name_not_in_list(std::vector<std::string> const& list_)
        : list(list_)
        { }
        
        template <class Matrix, class SymmGroup>
        bool operator() (Measurement_Term<Matrix, SymmGroup> const& mterm) const
        {
            return std::find(list.begin(), list.end(), mterm.name) == list.end();
        }
        
    private:
        std::vector<std::string> const& list;
    };
}

// implement me for your purpose!
template<class Matrix, class SymmGroup>
class Measurements
{
public:
    typedef Measurement_Term<Matrix, SymmGroup> mterm_t;
    typedef typename mterm_t::type_t meas_type_t;
    typedef typename mterm_t::op_t op_t;
    
    enum system_type_t {Wavefunction, Densitymatrix};
    
    Measurements() : super_meas(false) { }
    
    Measurements(std::vector<op_t> const& ident, std::vector<op_t> const& fill,
                 system_type_t system_type=Wavefunction)
    : super_meas( system_type == Densitymatrix )
    , identities_(ident)
    , fillings_(fill)
    { }
    
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
    
    Measurements sublist(std::vector<std::string> const& meas_list) const
    {
        Measurements sublist(*this);
        sublist.terms.erase_if( detail::name_not_in_list(meas_list) );
        return sublist;
    }
    
    const std::vector<op_t>& identity_matrices() const
    {
    	return identities_;
    }
    const op_t& identity_matrix(int type) const
    {
    	return identity_matrices()[type];
    }
    
    const std::vector<op_t>& filling_matrices() const
    {
    	return fillings_;
    }
    const op_t& filling_matrix(int type) const
    {
    	return filling_matrices()[type];
    }
    
    void add_term (mterm_t const & term)
    {
    	terms.push_back(term.clone());
    }

    void clear ()
    {
    	terms.clear();
        identities_.clear();
        fillings_.clear();
    }
    
    bool is_super_meas() const
    {
        return super_meas;
    }
    
protected:
    bool super_meas;
    std::vector<op_t> identities_, fillings_;
    boost::ptr_vector<mterm_t> terms;
};


template <class Matrix, class SymmGroup>
void parse_overlaps(BaseParameters const & parms, size_t sweep, Measurements<Matrix, SymmGroup> & meas)
{
    /* Syntax for MEASURE_OVERLAP:
     *  (1) MEASURE_OVERLAP[obsname] = "/path/to/ckp.h5"
     *  (2) MEASURE_OVERLAP[obsname(sweep)] = "/path/to/ckp.h5"
     *
     * `obsname` is the name that will be given in archive output.
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
	os << " - Fill operator:" << std::endl << term.fill_operator << std::endl;
	os << " - Operators:" << std::endl << term.operators;
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
