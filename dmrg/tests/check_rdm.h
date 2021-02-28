/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2021         Leon Freitag <lefreita@ethz.ch>
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
#ifndef CHECK_RDM_H
#define CHECK_RDM_H
#include "maquis_dmrg.h"

// Compare the (RDM) measurement element by element to the reference.
// Allow differences in sign (phase) in the whole matrix if allow_phase_difference==true
template<class Meas>
void check_measurement_mat(const Meas& meas, const Meas& reference, bool allow_phase_difference=false)
{
    double phase = 1.;
    for (auto&& it = reference.first.begin(); it != reference.first.end(); it++)
    {
        auto index = std::distance(meas.first.begin(), std::find(meas.first.begin(), meas.first.end(), *it));
        auto index_ref = std::distance(reference.first.begin(), it);
        // the phase can be 1 or -1, so set the phase in the first run
        if (it == reference.first.begin() && allow_phase_difference)
            phase = std::copysign(1., meas.second[index]/reference.second[index_ref]);

        BOOST_CHECK_SMALL(std::abs(meas.second[index]-phase*reference.second[index_ref]), 5e-7);
    }
}

#endif