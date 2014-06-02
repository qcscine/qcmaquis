/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Sebastian Keller <sebkelleb@phys.ethz.ch>
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

#ifndef SYMMETRY_NU1LPG_H
#define SYMMETRY_NU1LPG_H

#include "utils/io.hpp"
#include <vector>
#include <list>

#include <boost/lexical_cast.hpp>
#include <boost/functional/hash.hpp>

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/array.hpp>

#include <alps/numeric/matrix.hpp>

#include <dmrg/block_matrix/symmetry/nu1pg.h>

template<int N, class S>
class NU1LPG;

template<int N, class S = int>
class NU1LPG
{
public:
    typedef S subcharge;
    typedef NU1ChargePG<N, S> charge;
    typedef std::vector<charge> charge_v;
    
    static const charge IdentityCharge;
    static const bool finite = false;
    static const alps::numeric::matrix<S> mult_table;

    static charge fuse(charge a, charge b)
    {
        return plus<NU1LPG, N, S>(a,b);
    }
    
    template<int R> static charge fuse(boost::array<charge, R> const & v)
    {
        charge ret = v[0];
        for (int i = 1; i < R; ++i)
            ret = fuse(ret, v[i]);
        return ret;
    }
};

template<class S>
alps::numeric::matrix<S> generate_large_mult_table()
{
    // ************* TO DO ***************
    // check double group --> insert input param in the function?
    // go in the right case for the double group
    // initialize variables:
    // number of irreps, mult table, inverse and adjoints elements?,
    //
    // ***********************************

    // Cinfv double group mapped to C64
    // inverse and adjoint elements not implemented --> sebastian didn't
    
    int num_irreps = 128;
    if(num_irreps/2 % 2 == 1){
        throw std::logic_error("Number of boson and fermion irreps must be even\n");}
    int shift = num_irreps/2;
    int irrep = 1;
    int mj = 0;
    std::vector<S> mj2rep(num_irreps+1);
    alps::numeric::matrix<S> mult_table(num_irreps,num_irreps);
    mj2rep[shift+mj] = irrep;
    
    // populate mj2rep vector with boson and fermion irreps
    for(mj = 2; mj <= num_irreps/2-2; mj+=2){
        irrep++;
        mj2rep[shift+mj] = irrep;
        irrep++;
        mj2rep[shift-mj] = irrep;
    }

    mj = num_irreps/2;
    irrep++;
    mj2rep[shift+mj] = irrep;
    mj2rep[shift-mj] = irrep;

    for(mj = 1; mj <= num_irreps/2-1; mj+=2){
        irrep++;
        mj2rep[shift+mj] = irrep;
        irrep++;
        mj2rep[shift-mj] = irrep;
    }

    // build multiplication table
    int mij = 0;
    int jrrep = 0;
    int ijrrep = 0;
    for(int mi = -num_irreps/2; mi <= num_irreps/2; mi++){
        for(mj = -num_irreps/2; mj <= num_irreps/2; mj++){
            mij = mi + mj;
            if(mij <  -num_irreps/2){mij = mij + num_irreps;}
            if(mij >   num_irreps/2){mij = mij - num_irreps;}
            if(mij == -num_irreps/2){mij = num_irreps/2;}
            irrep  = mj2rep[shift+mi];
            jrrep  = mj2rep[shift+mj];
            ijrrep = mj2rep[shift+mij];
            mult_table(irrep-1,jrrep-1) = ijrrep;
        }
    }

    //for(int ii=0; ii < num_irreps; ++ii){
    //    for(int jj=0; jj < num_irreps; ++jj){
    //        std::cout << mult_table(ii,jj) << "\t";
    //    }
    //    std::cout << std::endl;
    //}

    return mult_table;
    /*
    //old pg matrix of sebastian
    alps::numeric::matrix<S> r(8,8);
    r(0,0) = 0; r(0,1) = 1; r(0,2) = 2; r(0,3) = 3;   r(0,4) = 4; r(0,5) = 5; r(0,6) = 6; r(0,7) = 7;
    r(1,0) = 1; r(1,1) = 0; r(1,2) = 3; r(1,3) = 2;   r(1,4) = 5; r(1,5) = 4; r(1,6) = 7; r(1,7) = 6;
    r(2,0) = 2; r(2,1) = 3; r(2,2) = 0; r(2,3) = 1;   r(2,4) = 6; r(2,5) = 7; r(2,6) = 4; r(2,7) = 5;
    r(3,0) = 3; r(3,1) = 2; r(3,2) = 1; r(3,3) = 0;   r(3,4) = 7; r(3,5) = 6; r(3,6) = 5; r(3,7) = 4;

    r(4,0) = 4; r(4,1) = 5; r(4,2) = 6; r(4,3) = 7;   r(4,4) = 0; r(4,5) = 1; r(4,6) = 2; r(4,7) = 3;
    r(5,0) = 5; r(5,1) = 4; r(5,2) = 7; r(5,3) = 6;   r(5,4) = 1; r(5,5) = 0; r(5,6) = 3; r(5,7) = 2;
    r(6,0) = 6; r(6,1) = 7; r(6,2) = 4; r(6,3) = 5;   r(6,4) = 2; r(6,5) = 3; r(6,6) = 0; r(6,7) = 1;
    r(7,0) = 7; r(7,1) = 6; r(7,2) = 5; r(7,3) = 4;   r(7,4) = 3; r(7,5) = 2; r(7,6) = 1; r(7,7) = 0;
    return r;
    */
}

template<int N, class S> const typename NU1LPG<N,S>::charge NU1LPG<N,S>::IdentityCharge = typename NU1PG<N,S>::charge();
template<int N, class S> const alps::numeric::matrix<S> NU1LPG<N,S>::mult_table = generate_large_mult_table<S>();

typedef NU1LPG<2> TwoU1LPG;

#endif
