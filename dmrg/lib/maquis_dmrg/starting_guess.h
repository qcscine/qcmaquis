/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
*               2019         Stefan Knecht <stknecht@ethz.ch>
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
#ifndef STARTING_GUESS_H
#define STARTING_GUESS_H

#include "maquis_dmrg.h"

namespace maquis
{

    // Interface that runs preliminary calculations for several states to use them for Fiedler ordering and/or CI-DEAS starting guess
    template <class V>
    class StartingGuess
    {
        public:

            // Constructor
            // parms: common parameters for all states
            //
            StartingGuess(const DmrgParameters& parms, int nstates, const std::string & pname, bool do_fiedler, bool do_cideas,
                          const std::vector<std::vector<int> >& hf_occupations = std::vector<std::vector<int> >());
            ~StartingGuess();

            // calculate Fiedler order
            std::string getFiedlerOrder();

            // Calculate and save MPS for CI-DEAS guess, as "pname.checkpoint_state.X.h5"
            // Currently we do not return the MPS but rather save them into files, which can be read later from disk
            void cideas();

        private:

            class Impl;
            std::unique_ptr<Impl> impl_;
    };
}

#endif
