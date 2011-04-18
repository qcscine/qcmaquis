#include "models.h"

#include <alps/parameters.h>
#include "alps_lattice.hpp"
#include "hamiltonians.hpp"
#include "alps_model.hpp"

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
    
    void model_factory (std::string const & type, std::ifstream & ifs,
                        Lattice & lattice, Hamiltonian & H)
    {
        
        if (type == "alps") {
            alps::Parameters parms(ifs);
            ALPSLattice lattice_tmp(parms);
            ALPSModel H_tmp(lattice_tmp, parms);
            lattice = lattice_tmp;
            H = H_tmp;
        } else {
            throw std::runtime_error("Don't know this type of lattice / model!");
        }
        
    }
    
    
    
} // namespace
