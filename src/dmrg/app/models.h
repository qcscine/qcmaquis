#ifndef APP_MODELS_H
#define APP_MODELS_H

#include "lattice.h"
#include "hamiltonian.h"
#include "measurements.h"
#include "utils/DmrgParameters.h"

namespace app {
    
    template <class Matrix, class SymmGroup>
    ModelParameters model_parser (std::string const & type, std::string const & fname,
                       Lattice* & lattice, Hamiltonian<Matrix, SymmGroup>* & H,
                       typename SymmGroup::charge& initc, Measurements<Matrix, SymmGroup>& meas);
}


#endif
