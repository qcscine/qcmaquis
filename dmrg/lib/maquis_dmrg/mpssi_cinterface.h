/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
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

    // get 1-TDM. Disabled for now since it isn't needed in MPSSI
    // void qcmaquis_mpssi_get_onetdm(char* bra_pname, int bra_state, char* ket_pname, int ket_state, int* indices, V* values, int size);

    // get spin-components of tdm
    // input: bra_pname and ket_pname: project names for bra and ket
    //        bra_state and ket_state: state indices
    // size: size of tdmaa and tdmbb matrices
    // output: tdmaa -- non-symmetric TDM with spin-up spin-up
    //         tdmbb -- idem with spin-down spin-down
    void qcmaquis_mpssi_get_onetdm_spin(char* bra_pname, int bra_state, char* ket_pname, int ket_state,  V* tdmaa, V* tdmbb, int size);
    void qcmaquis_mpssi_transform(char* pname, int state, int Ms);
    void qcmaquis_mpssi_rotate(char* pname, int state, V* t, int t_size, V scale_inactive, int Ms);

}

#endif