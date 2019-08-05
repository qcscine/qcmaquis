/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef DUMP_XVECTOR_H
#define DUMP_XVECTOR_H
#include "dmrg/mp_tensors/xvector.h"

namespace measurements {
    // A simple measurement that dumps a vector with nonredundant MPS parameters
    // for a given MPS
    // The variation should be zero for an MPS with itself as a reference, so we will be dumping
    // zeros initially, but we need them to obtain the total number of MPS parameters
    template <class Matrix, class SymmGroup>
    class DumpXVector : public measurement<Matrix, SymmGroup> {
        typedef measurement<Matrix, SymmGroup> base;
    public:

        DumpXVector(const std::string& name_, const std::string& chkpname, const std::string& aux_filename) : base(name_), chkpname_(chkpname), aux_filename_(aux_filename) {}

        void evaluate(MPS<Matrix, SymmGroup> const& mps, boost::optional<reduced_mps<Matrix, SymmGroup> const&> rmps = boost::none)
        {
            // we construct the XVector and dump it into a text file (for Molcas) and a checkpoint
            lr::XVector<Matrix, SymmGroup> x(mps, mps);
            x.save(chkpname_);
            x.dump_to_textfile(aux_filename_);
        }

    protected:
        measurement<Matrix, SymmGroup>* do_clone() const
        {
            return new DumpXVector(*this);
        }
    private:
        std::string chkpname_;
        std::string aux_filename_;
    };

}
#endif