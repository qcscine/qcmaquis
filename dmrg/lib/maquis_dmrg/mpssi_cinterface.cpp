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

#include "mpssi_cinterface.h"

std::unique_ptr<maquis::MPSSIInterface<double> > mpssi_interface_ptr;

extern "C"
{
    typedef double V;
    void qcmaquis_mpssi_init(char* project_names[], int* states, int* nstates, int n_projects)
    {
        std::vector<std::vector<int> > states_(n_projects);
        int nstate_counter = 0;

        for (int i = 0; i < states_.size(); i++)
        {
            states_[i].resize(nstates[i]);
            std::copy(states+nstate_counter, states+nstate_counter+nstates[i], states_[i].begin());
            nstate_counter += nstates[i];
        }

        std::vector<std::string> project_names_(project_names, project_names + n_projects);
        mpssi_interface_ptr.reset(new maquis::MPSSIInterface<double>(project_names_, states_));
    }

    V qcmaquis_mpssi_overlap(char* bra_pname, int bra_state, char* ket_pname, int ket_state, bool su2u1)
    {
        std::string bra_pname_(bra_pname);
        std::string ket_pname_(ket_pname);
        return mpssi_interface_ptr->overlap(bra_pname_, bra_state, ket_pname_, ket_state, su2u1);
    }

    /*
    void qcmaquis_mpssi_get_onetdm(char* bra_pname, int bra_state, char* ket_pname, int ket_state, int* indices, V* values, int size)
    {
        std::string bra_pname_(bra_pname);
        std::string ket_pname_(ket_pname);

        // copy-paste from 1-RDM code from maquis_cinterface.cpp
        const typename maquis::meas_with_results_type<V>& meas = mpssi_interface_ptr->onetdm(bra_pname_, bra_state, ket_pname_, ket_state);
        // the size attribute is pretty much useless if we allocate the output arrays outside of the interface
        // we'll just use it to check if the size matches the size of the measurement
        //
        // When we have point group symmetry, OpenMOLCAS will allocate an array larger than our number of measurements
        // So as long as our measurements will fit into the allocated array, we're good
        // I.e. the sizes of an array do not have to match exactly
        assert(size >= meas.first.size());
        assert(size >= meas.second.size());
        for (int i = 0; i < meas.first.size(); i++)
        {
            values[i] = meas.second[i];
            indices[2*i] = meas.first[i][0];
            indices[2*i+1] = meas.first[i][1];
        }
    } */

    void qcmaquis_mpssi_get_onetdm_spin(char* bra_pname, int bra_state, char* ket_pname, int ket_state, V* tdmaa, V* tdmbb, int size)
    {
        std::string bra_pname_(bra_pname);
        std::string ket_pname_(ket_pname);

        const std::vector<typename maquis::meas_with_results_type<V> >& meas = mpssi_interface_ptr->onetdm_spin(bra_pname_, bra_state, ket_pname_, ket_state);

        const int n_meas = 2; // only aa and bb are required, but can be changed easily to accomodate more measurements
        assert(meas.size() == n_meas); // two measurements for aa, bb, in this order
        bool bra_eq_ket = ((bra_pname_ == ket_pname_) && (bra_state == ket_state));

        assert(size >= meas[0].first.size());
        assert(size >= meas[0].second.size());

        // find the maximum index to deduce the number of orbitals
        std::vector<int> idx_vector(meas[0].first.size()); // first copy all the first indices into a separate vector
        std::transform(std::begin(meas[0].first), std::end(meas[0].first), idx_vector.begin(), [](const std::vector<int> & i){return i[0];});  // indices should be the same for all measurements, so we just take the ones from the first measurement (aa)
        int L = *std::max_element(idx_vector.cbegin(), idx_vector.cend())+1; // and then find the max index. Add 1 because orbital indices start from 0

        assert(L*L == size);

        // fill tdmaa
        for (int m = 0; m < meas[0].first.size(); m++)
        {
            int i = meas[0].first[m][0];
            int j = meas[0].first[m][1];
            tdmaa[L*i+j] = meas[0].second[m];
            if (bra_eq_ket) // symmetrise if bra==ket
                tdmaa[L*j+i] = tdmaa[L*i+j];
        }
        // and tdmbb
        for (int m = 0; m < meas[1].first.size(); m++)
        {
            int i = meas[1].first[m][0];
            int j = meas[1].first[m][1];
            tdmbb[L*i+j] = meas[1].second[m];
            if (bra_eq_ket)
                tdmbb[L*j+i] = tdmbb[L*i+j];
        }

        // Old code
        // copy the four measurement results into a single array
        // The layout of the array is as follows:
        // 0,0 aa 0,0 bb 0,1 aa .... 0,L bb 1,0 aa ... 1,L bb ... L,L bb
        // This is different from MOLCAS, which expects
        // 0,0 aa 0,0 ab 0,1 aa 0,1 ab ... 0,L aa 0,L ab 0,0 ba 0,0 bb ... 0,L ba 0,L bb 1,0 aa ... 1,L ab 1,0 ba ... 1,L bb ... L,L bb
        // but since we want to be compatible not only with MOLCAS, we will use the first layout, which is more generic
        // and store also the indexes
        // We expect the conversion to MOLCAS format at the MOLCAS side, despite the additional cost
        // (if you want to change this behaviour in the future, feel free to do so)
        // for (int m = 0; m < n_meas; m++) // loop through all measurements
        //     for (int s = 0; s < meas[m].first.size(); s++) // loop through all matrix elements
        //     {
        //         int i = meas[m].first[s][0];
        //         int j = meas[m].first[s][1];

        //         int C = n_meas*(L*i+j)+m;

        //         indices[n_meas*C] = i;
        //         indices[n_meas*C+1] = j;
        //         values[C] = meas[m].second[s];
        //     }
    }

    void qcmaquis_mpssi_transform(char* pname, int state)
    {
        std::string pname_(pname);
        maquis::transform(pname_, state);
    }

    void qcmaquis_mpssi_rotate(char* pname, int state, V* t, int t_size, V scale_inactive)
    {
        std::string pname_(pname);
        std::vector<V> t_vec(t, t + t_size);
        mpssi_interface_ptr->rotate(pname_, state, t_vec, scale_inactive);
    }

}