
#include "lattice.h"
#include "alps_lattice.hpp"


Lattice * lattice_factory (std::string const & lattice, std::ifstream & ifs)
{
    if (lattice == std::string("alps_lattice"))
        return new ALPSLattice(ifs);
    else {
        throw std::runtime_error("Don't know this lattice!");
        return NULL;
    }
}
