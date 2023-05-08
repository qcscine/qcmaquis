/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#include "dmrg/models/chem/2u1/model.h"
#include "dmrg/models/factories/factory.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1PG> {
    static std::shared_ptr<model_impl<Matrix, TwoU1PG> > parse
    (Lattice const & lattice, BaseParameters & parms)
    {
        typedef std::shared_ptr<model_impl<Matrix, TwoU1PG> > impl_ptr;
        if (parms["MODEL"] == std::string("quantum_chemistry")) {
            if (parms.is_set("LATTICE") && parms["LATTICE"] != std::string("orbitals"))
                throw std::runtime_error("Please use \"LATTICE = orbitals\" for quantum_chemistry\n");
            return impl_ptr( new qc_model<Matrix, TwoU1PG>(lattice, parms) );
        }
        else {
            throw std::runtime_error("Don't know this model: " + parms.get<std::string>("MODEL") + "\n");
            return impl_ptr();
        }
    }
};
