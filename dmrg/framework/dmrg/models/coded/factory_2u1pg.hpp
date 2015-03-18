/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#include "dmrg/models/chem/2u1/model_qc.h"

template<class Matrix>
struct coded_model_factory<Matrix, TwoU1PG> {
    static boost::shared_ptr<model_impl<Matrix, TwoU1PG> > parse
    (Lattice const & lattice, BaseParameters & parms)
    {
        typedef boost::shared_ptr<model_impl<Matrix, TwoU1PG> > impl_ptr;
        if (parms["MODEL"] == std::string("quantum_chemistry")) {
            if (parms["LATTICE"] != std::string("orbitals"))
                throw std::runtime_error("Please use \"LATTICE = orbitals\" for quantum_chemistry\n");

            return impl_ptr( new qc_model<Matrix, TwoU1PG>(lattice, parms) );
        }

        else {
            throw std::runtime_error("Don't know this model: " + parms.get<std::string>("MODEL") + "\n");
            return impl_ptr();
        }
    }
};
