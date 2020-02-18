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
    void qcmaquis_mpssi_init(int nel, int* multiplicities, int n_multiplicities,
                             int* states, int* nstates, char* pname, char* mult_suffixes[])
    {
        std::vector<int> multiplicities_(multiplicities, multiplicities + n_multiplicities);
        std::vector<std::vector<int> > states_(n_multiplicities);
        int nstate_counter = 0;

        for (int i = 0; i < states_.size(); i++)
        {
            states_[i].resize(nstates[i]);
            std::copy(states+nstate_counter, states+nstate_counter+nstates[i], states_[i].begin());
            nstate_counter += nstates[i];
        }

        std::string pname_(pname);
        std::vector<std::string> mult_suffixes_(mult_suffixes, mult_suffixes + n_multiplicities);
        mpssi_interface_ptr.reset(new maquis::MPSSIInterface<double>(nel, multiplicities_,
                                  states_, pname_, mult_suffixes_));
    }

    V qcmaquis_mpssi_overlap(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity, bool su2u1)
    {
        return mpssi_interface_ptr->overlap(bra_state, bra_multiplicity, ket_state, ket_multiplicity, su2u1);
    }

    void qcmaquis_mpssi_get_onetdm(int bra_state, int bra_multiplicity, int ket_state, int ket_multiplicity, int* indices, V* values, int size)
    {
        // copy-paste from 1-RDM code from maquis_cinterface.cpp
        const typename maquis::meas_with_results_type<V>& meas = mpssi_interface_ptr->onetdm(bra_state, bra_multiplicity, ket_state, ket_multiplicity);
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
    }

    void qcmaquis_mpssi_transform(char* checkpoint_name, char* suffix, int state, int multiplicity)
    {
        std::string pname_(checkpoint_name);
        std::string suffix_(suffix);
        mpssi_interface_ptr->transform(pname_, suffix, state, multiplicity);
    }
    void qcmaquis_mpssi_rotate(char* checkpoint_name, V* t, int t_size, V scale_inactive)
    {
        std::string pname_(checkpoint_name);
        std::vector<V> t_vec(t, t + t_size);
        mpssi_interface_ptr->rotate(pname_, t_vec, scale_inactive);
    }

}