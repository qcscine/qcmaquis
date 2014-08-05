/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef ENGINES_H
#define ENGINES_H

// Move this to separate compilation unit once done

#include <boost/shared_ptr.hpp>

#include "dmrg/mp_tensors/contractions.h"

template <class Matrix, class OtherMatrix, class SymmGroup>
boost::shared_ptr<contraction::Engine<Matrix, OtherMatrix, SymmGroup> >
engine_factory(BaseParameters & parms)
{
#ifdef ENABLE_SU2
    typedef boost::shared_ptr<contraction::Engine<Matrix, OtherMatrix, SymmGroup> > impl_ptr;
    if (parms["MODEL"] == "quantum_chemistry_SU2")
        return impl_ptr(new contraction::SU2Engine<Matrix, OtherMatrix, SymmGroup>());
    else if (parms["MODEL"] == "fermion Hubbard SU2")
        return impl_ptr(new contraction::SU2Engine<Matrix, OtherMatrix, SymmGroup>());
    else
        return impl_ptr(new contraction::AbelianEngine<Matrix, OtherMatrix, SymmGroup>());

#else
    return impl_ptr(new contraction::AbelianEngine<Matrix, OtherMatrix, SymmGroup>());

#endif
}

#endif
