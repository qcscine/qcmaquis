/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/MolecularHamiltonians/2u1/model.h"
#include "dmrg/models/factories/factory.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1PG> {
public:
    /** @brief Parses input info and returns a pointer to the proper model */
    static std::shared_ptr<model_impl<Matrix, TwoU1PG> > parse(const Lattice& lattice, BaseParameters& parms)
    {
        using impl_ptr = std::shared_ptr<model_impl<Matrix, TwoU1PG>>;
        if (parms["MODEL"] == std::string("quantum_chemistry")) {
            if (parms.is_set("LATTICE") && parms["LATTICE"] != std::string("orbitals"))
                throw std::runtime_error("Please use \"LATTICE = orbitals\" for quantum_chemistry\n");
            if (parms.is_set("transcorrelated_hamiltonian"))
                return (parms["transcorrelated_hamiltonian"] == true) ?
                        impl_ptr( new qc_model<Matrix, TwoU1PG, Hamiltonian::Electronic, HamiltonianTransformation::Transcorrelated>(lattice, parms) ) :
                        impl_ptr( new qc_model<Matrix, TwoU1PG, Hamiltonian::Electronic, HamiltonianTransformation::Conventional>(lattice, parms) );
            else
                return impl_ptr( new qc_model<Matrix, TwoU1PG, Hamiltonian::Electronic, HamiltonianTransformation::Conventional>(lattice, parms) );
        }
        else {
            throw std::runtime_error("Don't know this model: " + parms.get<std::string>("MODEL") + "\n");
            return impl_ptr();
        }
    }
};
