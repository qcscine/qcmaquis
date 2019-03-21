/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018-2019    Leon Freitag <lefreita@ethz.ch>
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
#ifndef MAQUIS_DMRG_DETAIL_H
#define MAQUIS_DMRG_DETAIL_H

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/sim/symmetry_factory.h"
#include "dmrg/sim/abstract_sim.h"

// #include <alps/numeric/matrix/matrix.hpp>
// #include <complex>
// typedef alps::numeric::matrix<double> matrix;
// typedef alps::numeric::matrix<std::complex<double> > cmatrix;

struct simulation_base {
    virtual ~simulation_base() {}
    virtual void run() =0;
    virtual void measure() =0;
};

template <class SymmGroup>
class simulation : public simulation_base {
    public:
        simulation(DmrgParameters & parms_);
        void run();
        void measure();
    private:
        DmrgParameters & parms;
        std::shared_ptr<abstract_interface_sim> sim_ptr;

};

struct simulation_traits {
    typedef std::shared_ptr<simulation_base> shared_ptr;
    template <class SymmGroup> struct F {
        typedef simulation<SymmGroup> type;
    };
};

#endif