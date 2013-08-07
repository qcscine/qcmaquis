/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_HAMILTONIANS_H
#define QC_HAMILTONIANS_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <boost/shared_ptr.hpp>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/term_maker.h"
#include "dmrg/models/chem/chem_detail.h"


template<class Matrix>
class qc_model : public Model<Matrix, TwoU1>, HamiltonianTraits
{
    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:
    
    qc_model(Lattice const & lat_, BaseParameters & parms_) : lat(lat_), parms(parms_)
    {
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
        
    }

    Index<TwoU1> get_phys() const
    {
        return phys;
    }
                            
    Hamiltonian<Matrix, TwoU1> H () const
    {
        return H_impl<Matrix>();
    }

    /* Disabled - need to implement iterators for one_matrix */
    /*
    Hamiltonian<one_matrix, TwoU1> H_chem () const
    {
        return H_impl<one_matrix>();
    }
    */
    
    Measurements<Matrix, TwoU1> measurements () const
    {
        return Measurements<Matrix, TwoU1>();
    }

    template <class M>
    Hamiltonian<M, TwoU1> H_impl () const;

private:
    Index<TwoU1> phys;

    Lattice const & lat;
    BaseParameters & parms;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }

};


#include "dmrg/models/chem/model_qc.hpp"

#endif
