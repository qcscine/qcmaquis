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
    // project_names: array of project names (array of size n_projects)
    // states: state indexes to be considered for each project
    // nstates: number of states considered for each project (array of the size n_projects)
    void qcmaquis_mpssi_init(char* project_names[], int* states, int* nstates, int n_projects);

    // Calculate overlap
    V qcmaquis_mpssi_overlap(char* bra_pname, int bra_state, char* ket_pname, int ket_state, bool su2u1);

    //
    void qcmaquis_mpssi_get_onetdm(char* bra_pname, int bra_state, char* ket_pname, int ket_state, int* indices, V* values, int size);

    // get spin-components of tdm
    // input: bra_pname and ket_pname: project names for bra and ket
    //        bra_state and ket_state: state indices
    // size: size of tdmaa and tdmbb matrices
    // output: tdmaa -- non-symmetric TDM with spin-up spin-up
    //         tdmbb -- idem with spin-down spin-down
    void qcmaquis_mpssi_get_onetdm_spin(char* bra_pname, int bra_state, char* ket_pname, int ket_state,  V* tdmaa, V* tdmbb, int size);
    void qcmaquis_mpssi_transform(char* checkpoint_name, int state);
    void qcmaquis_mpssi_rotate(char* checkpoint_name, V* t, int t_size, V scale_inactive);

}

#endif