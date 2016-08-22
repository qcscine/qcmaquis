/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Laboratory of Physical Chemistry, ETH Zurich
 *               2016 by Stefan Knecht <stknecht@ethz.ch>
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

#include <iostream>
#include <algorithm>
#include <string>
#include <boost/lexical_cast.hpp>

#include "dmrg/sim/matrix_types.h"
#include "dmrg/models/model.h"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/models/lattice.h"
#include "alps/numeric/matrix.hpp"
#include "dmrg/models/chem/util.h"
#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"
#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps_join.h"
#include "dmrg/utils/storage.h"

#include "dmrg/models/chem/transform_symmetry.hpp"
#include "../tools/ci_encode.hpp"
#include "utils.hpp"


#ifdef USE_AMBIENT
    #include "dmrg/block_matrix/detail/ambient.hpp"
    typedef ambient::numeric::tiles<ambient::numeric::matrix<double> > matrix;
#else
    #include "dmrg/block_matrix/detail/alps.hpp"
    typedef alps::numeric::matrix<double> matrix;
#endif


#if defined(USE_TWOU1)
typedef TwoU1 grp;
#elif defined(USE_TWOU1PG)
typedef TwoU1PG grp;
#elif defined(USE_SU2U1)
typedef SU2U1 grp;
#elif defined(USE_SU2U1PG)
typedef SU2U1PG grp;
#elif defined(USE_NONE)
typedef TrivialGroup grp;
#elif defined(USE_U1)
typedef U1 grp;
#endif

// function to write MPS to file
template<class Matrix, class SymmGroup>
void dump_MPS(MPS<Matrix, SymmGroup> & mps,
              DmrgParameters & parms,
              std::string mps_in_file,
              int file_id)
{
    maquis::cout << "norm of MPS to be dumped: " << norm(mps) << std::endl; 
    std::string mps_out_file = mps_in_file;
    if(file_id >= 0){
        std::size_t pos = mps_out_file.find(".h5");
        if (pos != mps_out_file.size())
            mps_out_file.erase(pos, 3);
        mps_out_file += "." + boost::lexical_cast<std::string>(file_id) + ".h5";
        save(mps_out_file, mps);
        if (boost::filesystem::exists(mps_out_file + "/props.h5"))
            boost::filesystem::remove(mps_out_file + "/props.h5");
        boost::filesystem::copy(mps_in_file + "/props.h5", mps_out_file + "/props.h5");

        storage::archive ar_out(mps_out_file + "/props.h5", "w");
        ar_out["/parameters"] << parms;
    }
    else save(mps_out_file, mps);
}

// function to calculate MPS' = MPO|MPS> aka (in CI terminology) calculating the sigma vector: sigma = H*C
template<class Matrix, class SymmGroup>
MPS<Matrix, SymmGroup> MPS_sigma_vector_product(MPS<Matrix, SymmGroup> const & mps,
                                                MPO<Matrix, SymmGroup> const & mpo,
                                                DmrgParameters & parms)
{   
    //mpo_times_mps_contractor_ss<Matrix, SymmGroup, storage::nop> sigma_vector_product(mps, mpo, parms);
    //sigma_vector_product.sweep();

    MPS<Matrix, SymmGroup> MPS_sigma; 
    //return (MPS_sigma = sigma_vector_product.get_current_mps());
    return mps;
}

// function to set a new MPO with elements defined by an input integral file
template<class Matrix, class SymmGroup>
MPO<Matrix, SymmGroup> setupMPO(std::string file, Lattice const & lattice, DmrgParameters & parms)
{
    typedef int pos_t;
    maquis::cout << "reading integrals for MPO from file: " << file << std::endl;
    parms.set("integral_file", file);

    Model<Matrix,SymmGroup> model = Model<Matrix,SymmGroup>(lattice, parms);
    MPO<Matrix, SymmGroup> tmpMPO = make_mpo(lattice, model);
    for (pos_t j = 0; j < lattice.size(); ++j){
        maquis::cout << "MPO for site " << j << std::endl;
        for (int b1 = 0; b1 < tmpMPO[j].row_dim(); ++b1)
        {
            for (int b2 = 0; b2 < tmpMPO[j].col_dim(); ++b2)
            {
                 if (tmpMPO[j].has(b1, b2)){
                     maquis::cout << tmpMPO[j].tag_number(b1,b2) << " ";
                 }
                 else maquis::cout << ". ";
            }
            maquis::cout << std::endl;
        }
    }
    return tmpMPO;
}

// function to read a scaling factor from an input file
template<class Matrix, class SymmGroup>
typename Matrix::value_type get_scaling(std::string file)
{
    std::ifstream scaling_file;
    scaling_file.open(file.c_str());
    double scale_coeff;

    scaling_file >> scale_coeff;
    typename Matrix::value_type factor = boost::lexical_cast<double>(scale_coeff);
    scaling_file.close();
    return factor;
}

// scale MPS tensor according to physical charges
template<class Matrix, class SymmGroup>
void scale_MPSTensor(MPSTensor<Matrix, SymmGroup> & mps,
                     typename Matrix::value_type tjj)
{
    typedef std::size_t size_t;
    typedef typename SymmGroup::charge charge;

    mps.make_left_paired();
    //maquis::cout << mps << std::endl; 

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
                // this would be BLAS detachable, i.e.
                // boost::numeric::bindings::blas::detail::scal(lsize, tjj_local, &(*m1[block].col(k).first), 1);
                for (size_t j = 0; j < lsize; j++){
                    m1[block](in_left_offset+j, k) *= tjj_local;
                }
        }
    }
}

template <class Matrix, class SymmGroup>
void rotate_mps(MPS<Matrix, SymmGroup> & mps, std::string scale_fac_file, std::string fcidump_file, DmrgParameters & parms,
                Lattice & lattice)
{
    typedef Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    pos_t L = mps.length();

    // step 1: scale MPS wrt rotations among inactive orbitals
    value_type alpha = get_scaling<Matrix, SymmGroup>(scale_fac_file + "." + boost::lexical_cast<std::string>(0));
    mps[0].multiply_by_scalar(alpha);

    MPS<Matrix, SymmGroup> mps_intermediate, mps_prime, mps_prime_prime;

    // step 2: scale MPS wrt rotations among active orbitals
    for (pos_t j = 0; j < L; ++j)
    {
        maquis::cout << "ROTATION of site "<< j << std::endl << "---------------- "<<      std::endl; 

        // scale the j-th MPS tensor wrt the occupation of the j-th orbital 
        value_type tjj = get_scaling<Matrix, SymmGroup>(scale_fac_file + "." + boost::lexical_cast<std::string>(j+1));
        scale_MPSTensor<Matrix, SymmGroup>(mps[j], tjj);
                                                                    debug::mps_print(mps[j], "\nScaled MPS at site " + boost::lexical_cast<std::string>(j));

        // get MPO
        MPO<Matrix, SymmGroup> myMPO = setupMPO<Matrix, SymmGroup>(fcidump_file + "." + boost::lexical_cast<std::string>(j+1), 
                                                                   lattice, parms);
        // |mps'> = H|mps> (first correction vector)
        mps_prime = MPS_sigma_vector_product<Matrix, SymmGroup>(mps, myMPO, parms);
        for (pos_t k = 0; k < L; ++k)
            mps_prime[k].replace_left_paired(mps_prime[k].data());

                                                                    debug::mps_print(mps_prime, "First correction MPS at site ");
                                                                    debug::mps_print_ci(mps_prime, parms, "dets.txt");
        for (pos_t k = 0; k < L; ++k)
            mps[k].data() += mps_prime[k].data();
        for (pos_t k = 0; k < L; ++k)
            mps[k].replace_left_paired(mps[k].data());
                                                                    debug::mps_print(mps, "Intermediate MPS at site ");
        // |mps''> = H|mps'> (second correction vector)
        mps_prime_prime = MPS_sigma_vector_product<Matrix, SymmGroup>(mps_prime, myMPO, parms);
        for (pos_t k = 0; k < L; ++k)
            mps_prime_prime[k].replace_left_paired(mps_prime_prime[k].data());

                                                                    debug::mps_print(mps_prime_prime, "Second correction MPS at site ");
        // set new MPS := mps + mps' + 1/2 mps''
        mps_prime_prime[0].multiply_by_scalar(0.5);
        for (pos_t k = 0; k < L; ++k)
            mps[k].data() += mps_prime_prime[k].data();

        for (pos_t k = 0; k < L; ++k)
            mps[k].replace_left_paired(mps[k].data());
                                                                    debug::mps_print(mps, "Final (for the current site to be rotated) MPS at site ");
    }
}

int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << "dmrg-input" << std::endl;
            return 1;
        }

        typedef Lattice::pos_t pos_t;

        DmrgOptions opt(argc, argv);
        DmrgParameters parms = opt.parms;

        // set non-Hermitian MPO (no complex conjugate elements for t_ij)
        parms.set("HasCConjugate", false);

        std::string rfile(parms.get<std::string>("resultfile"));
        std::string mps_in_file(parms.get<std::string>("chkpfile"));
       
        Lattice lattice(parms);

        MPS<matrix, grp> mps;
        load(mps_in_file, mps);

        mps.canonize(0);
        debug::mps_print_ci(mps, parms, "dets.txt");
        //maquis::cout << "norm of MPS: " << norm(mps) << std::endl; 
        //debug::mps_print(mps, "Original MPS at site ");

        pos_t L = mps.length();             assert(L = lattice.size());

        std::string scale_fac_file  = "tjj.tramps.orb";
        std::string fcidump_file    = "FCIDUMP.tramps.orb";

        //rotate_mps(mps, scale_fac_file, fcidump_file, parms, lattice);

        //maquis::cout << " FINAL DATA" << std::endl << " ----------" << std::endl; 
        //debug::mps_print(mps, " Rotated MPS at site ");
        //maquis::cout << "norm of final MPS: " << norm(mps) << std::endl; 

        //dump_MPS<matrix, grp>(mps, parms, mps_in_file, -1);

        MPS<matrix, grp> rot_mps(L);
        MPO<matrix, grp> my_mpo = setupMPO<matrix, grp>(fcidump_file + "." + boost::lexical_cast<std::string>(1), lattice, parms);
        for (int j = 0; j < L; ++j)
        {
            grp::charge delta = grp::IdentityCharge;
            MPSTensor<matrix, grp> direct_product = mpo_times_mps(my_mpo[j], mps[j], delta);
            rot_mps[j] = direct_product;
            //rot_mps[j].make_left_paired();
            //rot_mps[j].make_right_paired();
            debug::mps_print(rot_mps[j], "");
        }
        maquis::cout << "\nNorm: " << norm(rot_mps) << std::endl;

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
