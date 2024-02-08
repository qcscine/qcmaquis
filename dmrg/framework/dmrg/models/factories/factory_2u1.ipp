/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/MolecularHamiltonians/2u1/model.h"
#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/factories/factory.h"
#include "dmrg/models/FermiHubbardModels/2u1/RealSpaceFermiHubbardModel.h"
#include "dmrg/models/FermiHubbardModels/2u1/MomentumSpaceFermiHubbardModel.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1> 
{
    // Types definition
    using PtrType = std::shared_ptr<model_impl<Matrix, TwoU1> >;
    static PtrType parse(Lattice const & lattice, BaseParameters & parms)
    {
        using impl_ptr = std::shared_ptr<model_impl<Matrix, TwoU1> >;
        if (parms["MODEL"] == std::string("quantum_chemistry")) {
            if (parms.is_set("transcorrelated_hamiltonian")) {
                return (parms["transcorrelated_hamiltonian"] == 1) ?
                        impl_ptr( new qc_model<Matrix, TwoU1, Hamiltonian::Electronic, HamiltonianTransformation::Transcorrelated>(lattice, parms) ) :
                        impl_ptr( new qc_model<Matrix, TwoU1, Hamiltonian::Electronic, HamiltonianTransformation::Conventional>(lattice, parms) );
            }
            else {
                return impl_ptr( new qc_model<Matrix, TwoU1, Hamiltonian::Electronic, HamiltonianTransformation::Conventional>(lattice, parms) );
            }
        }
        else if (parms["MODEL"] == std::string("fermi_hubbard_real")) {
            return (parms["transcorrelated_hamiltonian"] == 1) ?
                    impl_ptr(new FermiHubbardRealTwoU1<Matrix>(lattice, parms, true)) :
                    impl_ptr(new FermiHubbardRealTwoU1<Matrix>(lattice, parms, false));
        } else if (parms["MODEL"] == std::string("fermi_hubbard_momentum")) {
            return (parms["transcorrelated_hamiltonian"] == 1) ?
                    impl_ptr(new FermiHubbardMomentumTwoU1<Matrix>(lattice, parms, true)) :
                    impl_ptr(new FermiHubbardMomentumTwoU1<Matrix>(lattice, parms, false));
#if defined(HAVE_NU1) && defined(DMRG_PREBO)
        } else if (parms["MODEL"] == std::string("PreBO")) {
            return impl_ptr( new PreBO<Matrix, 2>(lattice, parms) );
#endif
#if defined(HAVE_NU1) && defined(DMRG_VIBRATIONAL)
        } else if (parms["MODEL"] == std::string("nmode")) {
            return impl_ptr( new NMode<Matrix, 2>(lattice, parms, false) );
#endif
        } else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
