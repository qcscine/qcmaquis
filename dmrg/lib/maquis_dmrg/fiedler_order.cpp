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
#include "fiedler_order.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/sim/matrix_types.h"


namespace maquis
{
    template <class V>
    struct FiedlerOrder<V>::Impl
    {

        std::string calculate_FiedlerOrder()
        {

            typedef alps::numeric::matrix<V> Matrix;

            Matrix rdm(2,2);
            rdm(0,0) = 1; rdm(1,0) = 1; rdm(0,1) = 1; rdm(1,1) = -1;
            Matrix evecs(2,2);
            std::vector<V> eval(2);
            alps::numeric::syev(rdm,evecs,eval);

            return "1,2,3,4";
        }
        Impl() = default;
       ~Impl() = default;
    };

    template <class V>
    FiedlerOrder<V>::FiedlerOrder(int nstates,
                                  const std::string& pname):
                                  nstates_(nstates),
                                  pname_(pname),
                                  impl_()
    {

        // add stuff here
        for (std::size_t i = 0; i < nstates_; i++)
          {
                /*
                std::string chkp_name = pname_ + "." + mult_suffixes_[TwoS] + ".checkpoint_state." +
                    std::to_string(states[TwoS][i]) + ".h5";

                // Transform all checkpoints into the 2U1 representation
                transform(pname, mult_suffixes_[TwoS], states_[TwoS][i], multiplicities_[TwoS]);

                // Convert SU2 name to 2U1 name
                std::string twou1_chkp_name = twou1_name(pname_, mult_suffixes_[TwoS], states[TwoS][i], multiplicities_[TwoS]);

                // Rotate if we have more than 1 multiplicity
                if (multiplicities_.size() > 1)
                    rotate(twou1_chkp_name);
                */
          }
    }

    template <class V>
    FiedlerOrder<V>::~FiedlerOrder() = default;

    template <class V>
    std::string FiedlerOrder<V>::calculate_FiedlerOrder()
    {
        return impl_->calculate_FiedlerOrder();
    }

    template class FiedlerOrder<double>;

}

/*
template <class Matrix>
{
    std::array<V,2> I;
    std::array<V,1> s1;
    double omega = 1.0/nstates_;
    for (int i = 0; i < nstates_; i++)
    {
        //interface_ptr->optimize();
        //const typename maquis::DMRGInterface<V>::meas_with_results_type& meas = interface_ptr->get_chem_entropy();

        // accumulate an average mutual information I and single-orbital entropy s1

    };

    // calculate the fiedler ordering from averaged values
    std::string fiedlerorder = calculate_FiedlerOrder(&I);

    // set the new orbital order
    parms.set("orbital_order", fiedlerorder);

}
*/

