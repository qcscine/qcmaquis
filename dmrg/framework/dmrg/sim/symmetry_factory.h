/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/symmetry.h"

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <map>
#include <string>

#include "utils/io.hpp"

#include "dmrg/utils/guess_symmetry.h"

#include <utils/maquis_mpi.h>

namespace dmrg {
    
    template <class TR>
    typename TR::shared_ptr symmetry_factory(DmrgParameters & parms)
    {
        typedef typename TR::shared_ptr ptr_type;
        std::map<std::string, ptr_type> factory_map;
        
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "This binary contains symmetries: ";
#ifdef HAVE_NU1
        factory_map["nu1"] = ptr_type(new typename TR::template F<NU1>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "nu1 ";
#endif
#ifdef HAVE_TrivialGroup
        factory_map["none"] = ptr_type(new typename TR::template F<TrivialGroup>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "none ";
#endif
#ifdef HAVE_U1
        factory_map["u1"] = ptr_type(new typename TR::template F<U1>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "u1 ";
#endif
#ifdef HAVE_U1DG
        factory_map["u1dg"] = ptr_type(new typename TR::template F<U1DG>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "u1dg ";
#endif
#ifdef HAVE_TwoU1
        factory_map["2u1"] = ptr_type(new typename TR::template F<TwoU1>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "2u1 ";
#endif
#ifdef HAVE_TwoU1PG
        factory_map["2u1pg"] = ptr_type(new typename TR::template F<TwoU1PG>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "2u1pg ";
#endif
#ifdef HAVE_Ztwo
        factory_map["Z2"] = ptr_type(new typename TR::template F<Ztwo>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "Z2 ";
#endif
#ifdef HAVE_SU2U1
        factory_map["su2u1"] = ptr_type(new typename TR::template F<SU2U1>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "su2u1 ";
#endif
#ifdef HAVE_SU2U1PG
        factory_map["su2u1pg"] = ptr_type(new typename TR::template F<SU2U1PG>::type());
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << "su2u1pg ";
#endif
        if(maquis::mpi__->getGlobalRank() == 0)
            maquis::cout << std::endl;
        
        
        std::string symm_name;
        if (!parms.is_set("symmetry")) {
#ifdef HAVE_NU1
            symm_name = "nu1";
#else
            if (parms["model_library"] == "alps")
                symm_name = guess_alps_symmetry(parms);
#endif
        } else {
            symm_name = parms["symmetry"].str();
        }
        
        if (factory_map.find(symm_name) != factory_map.end())
            return factory_map[symm_name];
        else
            throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");

        parallel::sync();
        return ptr_type();
    }

}
