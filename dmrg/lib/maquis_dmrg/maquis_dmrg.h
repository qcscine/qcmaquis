/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MAQUIS_DMRG_H
#define MAQUIS_DMRG_H

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/sim/symmetry_factory.h"

#include <boost/shared_ptr.hpp>

namespace maquis
{
    struct simulation_base {
        virtual ~simulation_base() {}
        virtual void run(DmrgParameters & parms) =0;
    };

    template <class SymmGroup>
    struct dmrg_simulation : public simulation_base {
        void run(DmrgParameters & parms);
    };

    template <class SymmGroup>
    struct measure_simulation : public simulation_base {
        void run(DmrgParameters & parms);
    };

    struct dmrg_simulation_traits {
        typedef boost::shared_ptr<simulation_base> shared_ptr;
        template <class SymmGroup> struct F {
            typedef dmrg_simulation<SymmGroup> type;
        };
    };

    struct measure_simulation_traits {
        typedef boost::shared_ptr<simulation_base> shared_ptr;
        template <class SymmGroup> struct F {
            typedef measure_simulation<SymmGroup> type;
        };
    };

} // maquis
#endif
