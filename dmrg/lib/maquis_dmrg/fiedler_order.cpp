/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#include <algorithm>
#include "maquis_dmrg.h"
#include "dmrg/sim/matrix_types.h"


typedef double V;


template <class V, class Matrix>
class FiedlerOrder
{
    std::unique_ptr<maquis::DMRGInterface<V> > interface_ptr;
    DmrgParameters parms;
    private:
        int nstates_;
    public:
    // constructor
    FiedlerOrder(int nstates) : nstates_(nstates) { }
    // copy constructor 
    FiedlerOrder(const FiedlerOrder& Copy) {*this = Copy;}

    void get_FiedlerOrder();
    std::string calculate_FiedlerOrder(std::array<V,2> *I);

};

template <class Matrix> 
std::string FiedlerOrder<V,Matrix>::calculate_FiedlerOrder(std::array<V,2> *I)
{
    Matrix rdm(2,2);
    rdm(0,0) = 1; rdm(1,0) = 1; rdm(0,1) = 1; rdm(1,1) = -1;
    Matrix evecs(2,2);
    std::vector<V> eval(2);
    alps::numeric::syev(rdm,evecs,eval);
    return "1,2,3,4";
}

template <class Matrix>
void FiedlerOrder<V,Matrix>::get_FiedlerOrder()
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

