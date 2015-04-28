/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2015-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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

#ifndef QC_CHEM_CHECKS_H
#define QC_CHEM_CHECKS_H

namespace chem_detail {

    template <class SymmGroup, class = void>
    struct check_irrep
    {
        void operator()(BaseParameters & parms, storage::archive & ar) { }
    };
    template <class SymmGroup>
    struct check_irrep<SymmGroup, typename boost::enable_if<symm_traits::HasPG<SymmGroup> > ::type>
    {
        void operator()(BaseParameters & parms, storage::archive & ar)
        {
            std::string file; ar["/parameters/chkpfile"] >> file;

            int iref;
            ar["/parameters/irrep"] >> iref;
            if( iref != parms["irrep"])
                throw std::runtime_error("The existing checkpoint file " + file + " has the wrong irrep\n");
        }
    };

    template <class SymmGroup, class = void>
    struct check_su2_qns
    {
        void operator()(BaseParameters & parms, storage::archive & ar) { }
    };
    template <class SymmGroup>
    struct check_su2_qns<SymmGroup, typename boost::enable_if<symm_traits::HasSU2<SymmGroup> > ::type>
    {
        bool operator()(BaseParameters & parms, storage::archive & ar)
        {
            std::string file; ar["/parameters/chkpfile"] >> file;

            int spin;
            ar["/parameters/spin"] >> spin;
            if( spin != parms["spin"])
                throw std::runtime_error("The existing checkpoint file " + file + " has the wrong spin\n");

            int nelec;
            ar["/parameters/nelec"] >> nelec;
            if( nelec != parms["nelec"])
                throw std::runtime_error("The existing checkpoint file " + file + " has a wrong number of electrons\n");
        }
    };

    template <class SymmGroup, class = void>
    struct check_2u1_qns
    {
        void operator()(BaseParameters & parms, storage::archive & ar) { }
    };
    template <class SymmGroup>
    struct check_2u1_qns<SymmGroup, typename boost::enable_if<symm_traits::Has2U1<SymmGroup> > ::type>
    {
        bool operator()(BaseParameters & parms, storage::archive & ar)
        {
            std::string file; ar["/parameters/chkpfile"] >> file;

            int u1;
            ar["/parameters/u1_total_charge1"] >> u1;
            if( u1 != parms["u1_total_charge1"])
                throw std::runtime_error("The existing checkpoint file " + file + " has the wrong number of alpha-electrons\n");

            int u2;
            ar["/parameters/u1_total_charge2"] >> u2;
            if( u2 != parms["u1_total_charge2"])
                throw std::runtime_error("The existing checkpoint file " + file + " has a wrong number of beta-electrons\n");
        }
    };

    template <class SymmGroup, class = void>
    struct check_u1_qns
    {
        void operator()(BaseParameters & parms, storage::archive & ar) { }
    };
    template <class SymmGroup>
    struct check_u1_qns<SymmGroup, typename boost::enable_if<symm_traits::HasU1<SymmGroup> > ::type>
    {
        bool operator()(BaseParameters & parms, storage::archive & ar)
        {
            std::string file; ar["/parameters/chkpfile"] >> file;

            int u1;
            ar["/parameters/u1_total_charge"] >> u1;
            if( u1 != parms["u1_total_charge"])
                throw std::runtime_error("The existing checkpoint file " + file + " has the wrong number of electrons\n");
        }
    };

    template <class SymmGroup>
    struct RestoreCheck
    {
        void operator()(BaseParameters & parms, storage::archive & ar)
        {
            std::string file; ar["/parameters/chkpfile"] >> file;

            std::string sym;
            ar["/parameters/symmetry"] >> sym;
            if (sym != parms["symmetry"])
            {
                throw std::runtime_error("The existing checkpoint file " + file + " has wrong symmetry group " + sym + "\n");
            }

            check_irrep<SymmGroup>()(parms, ar);
            check_u1_qns<SymmGroup>()(parms, ar);
            check_2u1_qns<SymmGroup>()(parms, ar);
            check_su2_qns<SymmGroup>()(parms, ar);
        }
    };

}

#endif
