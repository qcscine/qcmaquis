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

#ifndef ABELIAN_ENGINE_FACTORY_H
#define ABELIAN_ENGINE_FACTORY_H

#include <boost/shared_ptr.hpp>

#include "dmrg/mp_tensors/contractions/abelian/engine.hpp"

namespace contraction {

    template <class Matrix, class OtherMatrix, class SymmGroup>
    class AbelianEngineFactory : public EngineFactory<Matrix, OtherMatrix, SymmGroup>
    {
        typedef boost::shared_ptr<contraction::Engine<Matrix, OtherMatrix, SymmGroup> > engine_ptr;

    public:
        virtual engine_ptr makeEngine() { return engine_ptr(new AbelianEngine<Matrix, OtherMatrix, SymmGroup>()); }
    };

} // namespace contraction

#endif
