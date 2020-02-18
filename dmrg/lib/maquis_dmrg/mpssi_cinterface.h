/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2020         Leon Freitag <lefreita@ethz.ch>
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
#ifndef MPSSI_CINTERFACE_H
#define MPSSI_CINTERFACE_H

#include "mpssi_interface.h"

extern "C"
{
    typedef double V;

    // Init the MPSSI interface
    // Parms:
    // nel -- number of electrons
    // multiplicities -- array of multiplicities of size n_multiplicities
    // states -- flattened array of arrays with states for each multiplicity
    // nstates -- array of numbers of states of each multiplicity
    // pname -- project name
    // mult_suffixes -- checkpoint suffixes for each multiplicity to be added after the project name
    // n_mult_suffixes -- size of mult_suffixes array
    void qcmaquis_mpssi_init(int nel, int* multiplicities, int n_multiplicities,
                             int* states, int* nstates, char* pname, char* mult_suffixes[]);
    V qcmaquis_mpssi_overlap(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity, bool su2u1);
    void qcmaquis_mpssi_get_onetdm(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity, int* indices, V* values, int size);
    void qcmaquis_mpssi_transform(char* checkpoint_name, char* suffix, int state, int multiplicity);
    void qcmaquis_mpssi_rotate(char* checkpoint_name, V* t, int t_size, V scale_inactive);

}

#endif