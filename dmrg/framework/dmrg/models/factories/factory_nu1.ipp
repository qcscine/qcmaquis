/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/vibrational/nu1/model.hpp"
#include "dmrg/models/MolecularHamiltonians/2u1/model.h"
#include "dmrg/models/factories/factory.h"

template<class Matrix, int N>
struct coded_model_factory<Matrix, NU1_template<N>> {
    static std::shared_ptr<model_impl<Matrix, NU1_template<N>> > parse(Lattice const& lattice, BaseParameters & parms)
    {
        using impl_ptr = std::shared_ptr<model_impl<Matrix, NU1_template<N>>>;
        if (parms["MODEL"] == std::string("PreBO")) {
#ifdef DMRG_PREBO
            return impl_ptr( new PreBO<Matrix, N>(lattice, parms) );
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else if (parms["MODEL"] == std::string("nmode")) {
#ifdef DMRG_VIBRATIONAL
            return impl_ptr( new NMode<Matrix, N>(lattice, parms, false) );
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};


template<class Matrix>
struct coded_model_factory<Matrix, NU1_template<2>> {
    static std::shared_ptr<model_impl<Matrix, NU1_template<2>> > parse(Lattice const& lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, NU1_template<2>> > impl_ptr;
        if (parms["MODEL"] == std::string("PreBO")) {
#ifdef DMRG_PREBO
            return impl_ptr( new PreBO<Matrix, 2>(lattice, parms) );
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
        else if (parms["MODEL"] == std::string("nmode")) {
#ifdef DMRG_VIBRATIONAL
            return impl_ptr( new NMode<Matrix, 2>(lattice, parms, false) );
#else
            throw std::runtime_error("Don't know this model!");
#endif
        }
#if defined(HAVE_TwoU1)
        else if (parms["MODEL"] == std::string("quantum_chemistry")) {
            return impl_ptr( new qc_model<Matrix, TwoU1, Hamiltonian::Electronic, HamiltonianTransformation::Conventional>(lattice, parms) );
        }
#endif
        else {
            throw std::runtime_error("Don't know this model!");
            return impl_ptr();
        }
    }
};
