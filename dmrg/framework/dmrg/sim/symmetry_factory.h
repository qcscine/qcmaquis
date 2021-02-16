/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2012 by Michele Dolfi <dolfim@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
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

#ifndef SYMMETRY_FACTORY_H
#define SYMMETRY_FACTORY_H

#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/block_matrix/symmetry.h"

#undef tolower
#undef toupper
#include <boost/tokenizer.hpp>
#include <map>
#include <string>

#include "utils/io.hpp"

#include "dmrg/utils/guess_symmetry.h"
#include "dmrg/utils/checks.h"

namespace dmrg {

    template <class TR, typename... Args>
    typename TR::shared_ptr symmetry_factory(DmrgParameters & parms, Args& ... args)
    {
        typedef typename TR::shared_ptr ptr_type;

        ptr_type ret;

        std::string symm_name;
        if (!parms.is_set("symmetry")) {
#ifdef HAVE_NU1
            symm_name = "nu1";
#else
            if (parms["model_library"] == "alps") // TODO: should be removed eventually
                symm_name = guess_alps_symmetry(parms);
            else
            {
                // obtain symmetry from checkpoint if exists
                if (parms.is_set("chkpfile"))
                {
                    parms["symmetry"] = maquis::checks::detail::get_symmetry(parms["chkpfile"]);
                    symm_name = parms["symmetry"].str();
                }
                else
                    throw std::runtime_error("Symmetry is not set in the parameters and cannot be autodetected.");
            }
#endif
        } else {
            symm_name = parms["symmetry"].str();
        }

        maquis::cout << "This binary contains symmetries: ";
#ifdef HAVE_NU1
        maquis::cout << "nu1 ";
#endif
#ifdef HAVE_TrivialGroup
        maquis::cout << "none ";
#endif
#ifdef HAVE_U1
        maquis::cout << "u1 ";
#endif
#ifdef HAVE_U1DG
        maquis::cout << "u1dg ";
#endif
#ifdef HAVE_TwoU1
        maquis::cout << "2u1 ";
#endif
#ifdef HAVE_TwoU1PG
        maquis::cout << "2u1pg ";
#endif
#ifdef HAVE_Ztwo
        maquis::cout << "Z2 ";
#endif
#ifdef HAVE_SU2U1
        maquis::cout << "su2u1 ";
#endif
#ifdef HAVE_SU2U1PG
        maquis::cout << "su2u1pg ";
#endif
        maquis::cout << std::endl;

        if (symm_name == "") // dummy if in case all ifdefs are false (unlikely)
        {}
#ifdef HAVE_NU1
        else if (symm_name == "nu1")
            ret.reset(new typename TR::template F<NU1>::type(args...));
#endif
#ifdef HAVE_TrivialGroup
        else if (symm_name == "none")
            ret.reset(new typename TR::template F<TrivialGroup>::type(args...));
#endif
#ifdef HAVE_U1
        else if (symm_name == "u1")
            ret.reset(new typename TR::template F<U1>::type(args...));
#endif
#ifdef HAVE_U1DG
        else if (symm_name == "u1dg")
            ret.reset(new typename TR::template F<U1DG>::type(args...));
#endif
#ifdef HAVE_TwoU1
        else if (symm_name == "2u1")
            ret.reset(new typename TR::template F<TwoU1>::type(args...));
#endif
#ifdef HAVE_TwoU1PG
        else if (symm_name == "2u1pg")
            ret.reset(new typename TR::template F<TwoU1PG>::type(args...));
#endif
#ifdef HAVE_Ztwo
        else if (symm_name == "Z2")
            ret.reset(new typename TR::template F<Ztwo>::type(args...));
#endif
#ifdef HAVE_SU2U1
        else if (symm_name == "su2u1")
            ret.reset(new typename TR::template F<SU2U1>::type(args...));
#endif
#ifdef HAVE_SU2U1PG
        else if (symm_name == "su2u1pg")
            ret.reset(new typename TR::template F<SU2U1PG>::type(args...));
#endif
        else
            throw std::runtime_error("Don't know this symmetry group. Please, check your compilation flags.");

        return ret;

        parallel::sync();
        return ptr_type();
    }

}

#endif