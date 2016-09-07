/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2016 Laboratory of Physical Chemistry, ETH Zurich
 *               2016 by Stefan Knecht <stknecht@ethz.ch>
 *               2016 by Sebastian Keller <sebkelle@phys.ethz.ch>
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
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_sectors.h"
#include "dmrg/mp_tensors/mpo_times_mps.hpp"
#include "dmrg/mp_tensors/mps_join.h"
#include "dmrg/mp_tensors/compression.h"
#include "dmrg/mp_tensors/mpo.h"

#include "dmrg/models/model.h"
#include "dmrg/models/lattice.h"

#include "dmrg/utils/DmrgOptions.h"
#include "dmrg/utils/DmrgParameters.h"

#include "dmrg/models/generate_mpo.hpp"
#include "dmrg/models/chem/transform_symmetry.hpp"
#include "dmrg/models/chem/2u1/chem_helper.h"

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

// function to write MPS to file
template<class Matrix, class SymmGroup>
void dump_MPS(MPS<Matrix, SymmGroup> & mps,
              //DmrgParameters & parms,
              std::string mps_in_file,
              int file_id)
{
    //maquis::cout << "norm of MPS to be dumped: " << norm(mps) << std::endl; 
    std::string mps_out_file = mps_in_file;
    if(file_id >= 0){
        /*
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
        */
    }
    else save(mps_out_file, mps);
}

// function to calculate MPS' = MPO|MPS> aka (in CI terminology) calculating the sigma vector: sigma = H*C
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
        typename SymmGroup::charge delta = grp::IdentityCharge;
        MPS<Matrix, grp> product(mps.size());
        for (int p = 0; p < mps.size(); ++p)
            product[p] =  mpo_times_mps(mpo_vec[i][p], mps[p], delta);

        clean_mps(product);
        assert(sensible(product));

        ret = (i==0) ? product : join(ret, product);
        assert(sensible(ret));

        //debug::mps_print(ret, "intra product ");
        //debug::mps_print_ci(ret, "dets.txt");
    }
    return ret;
}

// function to set a new MPO with elements defined by an input integral file
template<class Matrix, class SymmGroup>
std::vector<MPO<Matrix, SymmGroup> > setupMPO(std::string file, size_t L, size_t Nup, size_t Ndown, std::string site_types)
{
    typedef Lattice::pos_t pos_t;
    typedef typename MPOTensor<Matrix, SymmGroup>::tag_type tag_type;
    typedef typename SymmGroup::subcharge sc_t;
    //maquis::cout << "reading integrals for MPO from file: " << file << std::endl;

    BaseParameters parms = chem_detail::set_2u1_parameters(L, Nup, Ndown);
    parms.set("integral_file", file);
    parms.set("integral_cutoff", 0.);
    parms.set("site_types", site_types);

    Lattice lat(parms);
    Model<Matrix,SymmGroup> model = Model<Matrix, SymmGroup>(lat, parms);

    std::vector<tag_type> ident, fill;
    for (size_t i = 0; i <= lat.maximum_vertex_type(); ++i)
    {
        ident.push_back(model.identity_matrix_tag(i));
        fill.push_back(model.filling_matrix_tag(i));
    }

    //chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, ident, fill, model.operators_table());
    //std::vector<typename Matrix::value_type> & matrix_elements = term_assistant.getMatrixElements(); 
    std::vector<typename Matrix::value_type> matrix_elements;
    alps::numeric::matrix<pos_t> idx;
    boost::tie(idx, matrix_elements) = chem_detail::parse_integrals<typename Matrix::value_type, SymmGroup>(parms, lat, false);

    std::vector<MPO<Matrix, SymmGroup> > ret;
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        std::vector<pos_t> positions;
        std::vector<tag_type> operators_up, operators_down;

        int i = idx(m, 0);
        int j = idx(m, 1);
        int k = idx(m, 2);
        int l = idx(m, 3);

        //maquis::cout << "integral term ("<< i <<"|" <<j<<") = "<< matrix_elements[m] << std::endl;

        assert( k==-1 && l==-1);
        positions.push_back(i);
        positions.push_back(j);
        operators_up.push_back(model.get_operator_tag("create_up", lat.get_prop<sc_t>("type", i)));
        operators_up.push_back(model.get_operator_tag("destroy_up", lat.get_prop<sc_t>("type", j)));
        operators_down.push_back(model.get_operator_tag("create_down", lat.get_prop<sc_t>("type", i)));
        operators_down.push_back(model.get_operator_tag("destroy_down", lat.get_prop<sc_t>("type", j)));

        ret.push_back(generate_mpo::make_1D_mpo(positions, operators_up, ident, fill, model.operators_table(), lat, matrix_elements[m]));
        ret.push_back(generate_mpo::make_1D_mpo(positions, operators_down, ident, fill, model.operators_table(), lat, matrix_elements[m]));

        /*
        MPO<Matrix, SymmGroup> const & mpo = *ret.rbegin();
        for (pos_t j = 0; j < lat.size(); ++j){
            maquis::cout << "MPOTensor for site " << j << std::endl;
            for (int b1 = 0; b1 < mpo[j].row_dim(); ++b1)
            {
                for (int b2 = 0; b2 < mpo[j].col_dim(); ++b2)
                {
                     if (mpo[j].has(b1, b2)){
                         maquis::cout << mpo[j].tag_number(b1,b2) << " ";
                     }
                     else maquis::cout << ". ";
                }
                maquis::cout << std::endl;
            }
        }
        */
    }

    return ret;
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

template <class Matrix, class SymmGroup>
void compress_mps(MPS<Matrix, SymmGroup> & mps, std::string text)
{
    maquis::cout << "- MPS compression - input MPS: "<< text << std::endl;

    matrix::value_type final_norm        = norm(mps);
    matrix::value_type compression_trace = 1.0;

    mps = compression::l2r_compress(mps, 8000, 1e-8, compression_trace);
    maquis::cout << "- compression trace          : "<< compression_trace << std::endl;
    mps[0].multiply_by_scalar(compression_trace*sqrt(final_norm));

}

template <class Matrix, class SymmGroup>
void rotate_mps(MPS<Matrix, SymmGroup> & mps, std::string scale_fac_file, std::string fcidump_file)
{
    typedef Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;

    typename SymmGroup::subcharge Ndown, Nup;
    pos_t L = mps.length();
    Nup = mps[L-1].col_dim()[0].first[0];
    Ndown = mps[L-1].col_dim()[0].first[1];
    std::string site_types = chem_detail::infer_site_types(mps);

    //maquis::cout << "- input MPS - "<<      std::endl;
    //debug::mps_print_ci(mps, "dets.txt");

    // step 1: scale MPS wrt rotations among inactive orbitals
    value_type alpha = get_scaling<Matrix, SymmGroup>(scale_fac_file + "." + boost::lexical_cast<std::string>(0));
    mps[0].multiply_by_scalar(alpha);

    //maquis::cout << "- scaled MPS (inactive orbitals) - "<<      std::endl;
    //debug::mps_print_ci(mps, "dets.txt");

    MPS<Matrix, SymmGroup> mps_prime, mps_prime_prime;

    // step 2: scale MPS wrt rotations among active orbitals
    for (pos_t j = 0; j < L; ++j)
    {
        maquis::cout << "ROTATION of site "<< j << std::endl << "---------------- "<<      std::endl;

        // scale the j-th MPS tensor wrt the occupation of the j-th orbital 
        value_type tjj = get_scaling<Matrix, SymmGroup>(scale_fac_file + "." + boost::lexical_cast<std::string>(j+1));
        scale_MPSTensor<Matrix, SymmGroup>(mps[j], tjj);

        //maquis::cout << "- scaled MPS (local sites) - "<<      std::endl;
        //debug::mps_print_ci(mps, "dets.txt");
        //debug::mps_print(mps[j], "\nScaled MPS at site " + boost::lexical_cast<std::string>(j));

        // get MPO
        std::vector<MPO<Matrix, SymmGroup> > MPO_vec
            = setupMPO<Matrix, SymmGroup>(fcidump_file + "." + boost::lexical_cast<std::string>(j+1), L, Nup, Ndown, site_types);
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

int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << "mps.h5" << std::endl;
            return 1;
        }

        typedef Lattice::pos_t pos_t;

        //std::string rfile(parms.get<std::string>("resultfile"));
       
        MPS<matrix, grp> mps;
        load(argv[1], mps);

        mps.canonize(0);

        //maquis::cout << " Input MPS: " << std::endl; 
        //debug::mps_print_ci(mps, "dets.txt");

        //maquis::cout << "norm of MPS: " << norm(mps) << std::endl; 
        //debug::mps_print(mps, "Original MPS at site ");

        std::string scale_fac_file  = "tjj.tramps.orb";
        std::string fcidump_file    = "FCIDUMP.tramps.orb";

        rotate_mps(mps, scale_fac_file, fcidump_file);

        matrix::value_type final_norm = norm(mps);
        //maquis::cout << "norm of final MPS: " << norm(mps) << std::endl; 

        // NOTE: this normalizes the final MPS and may invert signs
        //mps = compression::l2r_compress(mps, 10000, 1e-7);
        //mps[0].multiply_by_scalar(sqrt(final_norm));

        maquis::cout << " FINAL DATA" << std::endl << " ----------" << std::endl;
        //debug::mps_print(mps, " Rotated MPS at site ");
        //debug::mps_print_ci(mps, "dets.txt");

        maquis::cout << "norm of final MPS: " << norm(mps) << std::endl; 

        //dump_MPS<matrix, grp>(mps, parms, mps_in_file, -1);
        dump_MPS<matrix, grp>(mps, argv[1], -1);

    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
