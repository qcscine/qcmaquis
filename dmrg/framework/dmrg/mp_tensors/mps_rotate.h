/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Laboratory of Physical Chemistry, ETH Zurich
 *               2016 by Stefan Knecht <stknecht@ethz.ch>
 *               2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2019 by Leon Freitag <lefreita@ethz.ch>
 *               2021 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MPS_ROTATE_H
#define MPS_ROTATE_H

#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_sectors.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps_join.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/mpo.h"
#include "integral_interface.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/chem/transform_symmetry.hpp"
#include "dmrg/models/chem/su2u1/term_maker.h"

// Functions required for MPSSI
// Equation numbers are from S. Knecht et al, JCTC 2016, 12, 5881
namespace mps_rotate
{
    template<class Matrix, class SymmGroup>
    bool sensible(MPS<Matrix, SymmGroup> const & mps)
    {
        for (size_t p = 0; p < mps.size()-1; ++p)
            if (mps[p].col_dim() == mps[p+1].row_dim())
                continue;
            else
                return false;

        return true;
    }

    // scale MPS tensor according to physical charges
    // see eq. 44
    template<class Matrix, class SymmGroup>
    void scale_MPSTensor(MPSTensor<Matrix, SymmGroup> & mps,
                        typename Matrix::value_type tjj)
    {
        typedef std::size_t size_t;
        typedef typename SymmGroup::charge charge;

        mps.make_left_paired();
        maquis::cout << "scaling factor " << tjj << std::endl;

        Index<SymmGroup> const & physical_i = mps.site_dim();
        Index<SymmGroup> const & left_i     = mps.row_dim();

        block_matrix<Matrix, SymmGroup> & m1 = mps.data();
        ProductBasis<SymmGroup> in_left_pb(physical_i, left_i);

        // loop over blocks in MPS tensor (remember MPS tensor has a block matrix structure)
        for (size_t block = 0; block < m1.n_blocks(); ++block)
        {
            // loop over physical charges (aka local states); nothing needs to be done for charge <0,0>
            for (size_t i = 0; i < physical_i.size(); i++)
            {
                charge phys_charge = physical_i[i].first;
                if (phys_charge == SymmGroup::IdentityCharge) continue;

                size_t l = left_i.position(SymmGroup::fuse(m1.basis().left_charge(block), -phys_charge));
                // if left charge does not exist (meaning that there is no matching input charge in this block matrix)
                if(l == left_i.size()) continue;

                // provide the offset in the block matrix
                std::size_t in_left_offset = in_left_pb(phys_charge, left_i[l].first);
                std::size_t lsize = left_i[l].second;

                typename Matrix::value_type tjj_local = tjj;
                if(SymmGroup::particleNumber(phys_charge) == 2) tjj_local *= tjj;

                // scale each coefficient in the block-th matrix
                for (size_t k = 0; k < m1[block].num_cols(); k++)
                    // this would be BLAS dispatchable, i.e.
                    // boost::numeric::bindings::blas::detail::scal(lsize, tjj_local, &(*m1[block].col(k).first), 1);
                    for (size_t j = 0; j < lsize; j++){
                        m1[block](in_left_offset+j, k) *= tjj_local;
                    }
            }
        }
    }

    // function to set a new MPO with elements defined by integrals in parms
    // implemented as a functor for
    // setup MPO with elements as in Eq. 45
    // t: full transformation matrix
    // j: orbital index
    // lat, model : lattice and model for the operator
    template<class Matrix, class SymmGroup, class = void>
    struct setupMPO
    {
        std::vector<MPO<Matrix, SymmGroup> > operator()(const Matrix & t, int j, const Lattice& lat, const Model<Matrix, SymmGroup> & model)
        {
            typedef Lattice::pos_t pos_t;
            typedef typename MPOTensor<Matrix, SymmGroup>::tag_type tag_type;
            typedef typename SymmGroup::subcharge sc_t;
            std::vector<tag_type> ident, fill;
            for (int iSite = 0; iSite < lat.getMaxType(); iSite++)
            {
                ident.push_back(model.identity_matrix_tag(iSite));
                fill.push_back(model.filling_matrix_tag(iSite));
            }

            std::vector<MPO<Matrix, SymmGroup> > ret;
            for (int i=0; i < t.num_rows(); ++i)
                if (i != j)
                {
                        std::vector<pos_t> positions {i,j};
                        std::vector<tag_type> operators_up, operators_down;
                        operators_up.push_back(model.get_operator_tag("create_up", lat.get_prop<sc_t>("type", i)));
                        operators_up.push_back(model.get_operator_tag("destroy_up", lat.get_prop<sc_t>("type", j)));
                        operators_down.push_back(model.get_operator_tag("create_down", lat.get_prop<sc_t>("type", i)));
                        operators_down.push_back(model.get_operator_tag("destroy_down", lat.get_prop<sc_t>("type", j)));

                        ret.push_back(generate_mpo::make_1D_mpo(positions, operators_up, ident, fill, model.operators_table(), lat, t(i,j)/t(j,j)));
                        ret.push_back(generate_mpo::make_1D_mpo(positions, operators_down, ident, fill, model.operators_table(), lat, t(i,j)/t(j,j)));

                }

            return ret;
        }
    };

    /*
    // does not work
    template<class Matrix, class SymmGroup>
    struct setupMPO<Matrix, SymmGroup, symm_traits::enable_if_su2_t<SymmGroup>>
    {
        std::vector<MPO<Matrix, SymmGroup> > operator()(const Matrix & t, std::size_t j, const Lattice& lat, const Model<Matrix, SymmGroup> & model)
        {

            std::vector<MPO<Matrix, SymmGroup> > ret;

            // find the highest irreducible representation number
            // used to generate ops for all irreps 0..max_irrep
            typename SymmGroup::subcharge max_irrep = 0;
            for (Lattice::pos_t p=0; p < lat.size(); ++p)
                max_irrep = (lat.get_prop<typename SymmGroup::subcharge>("type", p) > max_irrep)
                                ? lat.get_prop<typename SymmGroup::subcharge>("type", p) : max_irrep;

            typename TermMakerSU2<Matrix, SymmGroup>::Operators ops = TermMakerSU2<Matrix, SymmGroup>::construct_operators(max_irrep, model.operators_table());
            typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection op_collection = TermMakerSU2<Matrix, SymmGroup>::construct_operator_collection(ops, max_irrep);

            for (std::size_t i = 0; i < t.num_rows(); ++i)
                if (i != j)
                {
                    std::vector<typename Model<Matrix,SymmGroup>::term_descriptor> terms;
                    // copy-paste from 1-RDM expval from tagged_nrankrdm.h
                    // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
                    terms.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                        true, op_collection.ident.no_couple, t(i,j)/t(j,j)*std::sqrt(2.), i, j, op_collection.create.couple_down, op_collection.create.fill_couple_up,
                                                            op_collection.destroy.couple_down, op_collection.destroy.fill_couple_up, lat
                    ));


                    // check if term is allowed by symmetry
                    // if(not measurements_details::checkpg<SymmGroup>()(terms[0], model.operators_table(), lat))
                    //         continue;

                    generate_mpo::TaggedMPOMaker<Matrix, SymmGroup> mpo_m(lat, op_collection.ident.no_couple, op_collection.ident_full.no_couple,
                                                                            op_collection.fill.no_couple, model.operators_table(), terms);
                    ret.push_back(mpo_m.create_mpo());
                }

                return ret;
        }
    };
    */

    // function to calculate MPS' = MPO|MPS> aka (in CI terminology) calculating the sigma vector: sigma = H*C
    // as used e.g. in Eqs. 46 and 48
    template<class Matrix, class SymmGroup>
    MPS<Matrix, SymmGroup> MPS_sigma_vector_product(MPS<Matrix, SymmGroup> const & mps,
                                                    std::vector<MPO<Matrix, SymmGroup> > const & mpo_vec)
    {

        //mpo_times_mps_contractor_ss<Matrix, SymmGroup, storage::nop> sigma_vector_product(mps, mpo, parms);
        //sigma_vector_product.sweep();
        assert(sensible(mps));
        MPS<Matrix, SymmGroup> ret;
        for (size_t i = 0; i < mpo_vec.size(); ++i)
        {
            //debug::mps_print(mps, "before round " + boost::lexical_cast<std::string>(i));
            typename SymmGroup::charge delta = SymmGroup::IdentityCharge;
            MPS<Matrix, SymmGroup> product(mps.size());
            for (int p = 0; p < mps.size(); ++p)
                product[p] =  MPOTimesMPSTraitClass<Matrix, SymmGroup>::mpo_times_mps_singleop(mpo_vec[i][p], mps[p], delta);

            clean_mps(product);
            assert(sensible(product));

            ret = (i==0) ? product : join(ret, product);
            assert(sensible(ret));

            //debug::mps_print(ret, "intra product ");
            //debug::mps_print_ci(ret, "dets.txt");
        }
        return ret;
    }

    // MPS compression to keep the dimensions reasonable
    template <class Matrix, class SymmGroup>
    void compress_mps(MPS<Matrix, SymmGroup> & mps, std::string text="")
    {
        maquis::cout << "- MPS compression - input MPS: "<< text << std::endl;

        typename Matrix::value_type final_norm        = norm(mps);
        typename Matrix::value_type compression_trace = 1.0;

        mps = compression::l2r_compress(mps, 8000, 1e-8, compression_trace);
        maquis::cout << "- compression trace          : "<< compression_trace << std::endl;
        mps[0].multiply_by_scalar(compression_trace*sqrt(final_norm));

    }

    // MPS rotation as described in Sections III.b.2.b and III.b.2.c
    // t: rotation matrix from Eqs. 35f but only for active orbitals
    // inactive_scaling: inactive scaling factor from Eq. 43
    template <class Matrix, class SymmGroup>
    void rotate_mps(MPS<Matrix, SymmGroup> & mps, const Matrix& t, typename Matrix::value_type inactive_scaling)
    {
        typedef Lattice::pos_t pos_t;
        typedef typename Matrix::value_type value_type;

        typename SymmGroup::subcharge Ndown, Nup;

        pos_t L = mps.length();

        assert(t.num_rows() == L);
        assert(t.num_cols() == L);

        Nup = mps[L-1].col_dim()[0].first[0];
        Ndown = mps[L-1].col_dim()[0].first[1];

        // Create a lattice and model to set up MPOs

        // create a fake integral map to build a model so that it's happy.
        // we don't use the Hamiltonian of the models so we don't care for these integrals,
        // but the Model constructor will fail if no integrals are defined
        const chem::integral_map<typename Matrix::value_type> fake_integrals = { { { 1, 1, 1, 1 },   0.0 } };

        BaseParameters parms = chem::detail::set_2u1_parameters(L, Nup, Ndown);
        parms.set("integrals_binary", chem::serialize(fake_integrals));
        parms.set("integral_cutoff", 0.);
        parms.set("site_types", chem::detail::infer_site_types(mps));
        Lattice lat(parms);
        Model<Matrix,SymmGroup> model = Model<Matrix, SymmGroup>(lat, parms);

        //maquis::cout << "- input MPS - "<<      std::endl;
        //debug::mps_print_ci(mps, "dets.txt");

        // step 1: scale MPS wrt rotations among inactive orbitals (Eq. 42)

        mps[0].multiply_by_scalar(inactive_scaling);

        //maquis::cout << "- scaled MPS (inactive orbitals) - "<<      std::endl;
        //debug::mps_print_ci(mps, "dets.txt");

        // step 2: scale MPS wrt rotations among active orbitals
        MPS<Matrix, SymmGroup> mps_prime, mps_prime_prime;
        for (pos_t j = 0; j < L; ++j)
        {
            maquis::cout << "ROTATION of site "<< j << std::endl << "---------------- "<<      std::endl;

            // scale the j-th MPS tensor wrt the occupation of the j-th orbital

            scale_MPSTensor<Matrix, SymmGroup>(mps[j], t(j,j));

            //maquis::cout << "- scaled MPS (local sites) - "<<      std::endl;
            //debug::mps_print_ci(mps, "dets.txt");
            //debug::mps_print(mps[j], "\nScaled MPS at site " + boost::lexical_cast<std::string>(j));

            // get MPO
            std::vector<MPO<Matrix, SymmGroup> > MPO_vec = setupMPO<Matrix, SymmGroup>()(t, j, lat, model);

            // check for non-zero MPO vector
            if(MPO_vec.size() == 0) continue;

            // |mps'> = H|mps> (first correction vector)
            mps_prime = MPS_sigma_vector_product<Matrix, SymmGroup>(mps, MPO_vec);

            maquis::cout << "- first correction MPS obtained - "<<      std::endl;
            //debug::mps_print_ci(mps_prime, "dets.txt");


            mps = join(mps, mps_prime);
            //debug::mps_print(mps, "Intermediate MPS at site ");
            maquis::cout << "- intermediate correction MPS obtained - "<<      std::endl;
            //debug::mps_print_ci(mps, "dets.txt");

            // compression of MPS'
            compress_mps<Matrix, SymmGroup>(mps_prime, "MPS prime");

            //maquis::cout << "- enter for second correction - "<<      std::endl;
            // |mps''> = H|mps'> (second correction vector)
            mps_prime_prime = MPS_sigma_vector_product<Matrix, SymmGroup>(mps_prime, MPO_vec);
            maquis::cout << "- Second correction MPS obtained - "<<      std::endl;
            //debug::mps_print_ci(mps_prime_prime, "dets.txt");
            //debug::mps_print(mps_prime_prime, "Second correction MPS at site ");

            // compression of MPS''
            //compress_mps<Matrix, SymmGroup>(mps_prime_prime, "MPS prime prime");

            // set new MPS := mps + mps' + 1/2 mps''
            mps_prime_prime[0].multiply_by_scalar(0.5);
            mps = join(mps, mps_prime_prime);
            //
            maquis::cout << "-  Final (for the current site to be rotated) MPS with no compression   - "<<      std::endl;
            //debug::mps_print_ci(mps, "dets.txt");

            // compression of final MPS
            compress_mps<Matrix, SymmGroup>(mps, "MPS final");

            maquis::cout << "-  Final (for the current site to be rotated) MPS with full compression - "<<      std::endl;
            //debug::mps_print_ci(mps, "dets.txt");
            //debug::mps_print(mps, "Final (for the current site to be rotated) MPS at site ");
        }
    }

}
#endif