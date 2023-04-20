/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
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
