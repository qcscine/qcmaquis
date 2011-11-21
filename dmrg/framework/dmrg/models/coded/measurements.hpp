
#ifndef MEASUREMENTS_HPP_
#define MEASUREMENTS_HPP_

#include "dmrg/models/measurements.h"

#include "dmrg/utils/BaseParameters.h"

namespace app {
    
    template<class Matrix, class SymmGroup>
    struct pre_measurements {
        void operator() (const Lattice& lattice, BaseParameters & model,
                         std::vector<Measurement_Term<Matrix, SymmGroup> >& terms,
                         typename Measurement_Term<Matrix, SymmGroup>::op_t& ident)
        { }
    };
    
    // include specializations of pre_measurements
#include "meas_u1.hpp"
#include "meas_2u1.hpp"
    
    template<class Matrix, class SymmGroup>
    class CodedMeasurements : public Measurements<Matrix, SymmGroup>
    {
        typedef Measurements<Matrix, SymmGroup> super_t;
    public:
        CodedMeasurements (const Lattice& lattice, BaseParameters & model)
        {
            pre_measurements<Matrix, SymmGroup>()(lattice, model, super_t::terms, super_t::ident);
        }
    };
    
    
} // namespace
#endif /* MEASUREMENTS_HPP_ */
