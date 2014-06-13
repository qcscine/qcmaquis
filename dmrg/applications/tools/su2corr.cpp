/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#include <cmath>
#include <iostream>
#include <iomanip>
#include <sys/time.h>
#include <sys/stat.h>

using std::cerr;
using std::cout;
using std::endl;

#include "dmrg/mp_tensors/contractions/su2.hpp"
#include "dmrg/models/chem/prepare_corr.hpp"

typedef alps::numeric::matrix<double> matrix;

typedef TwoU1PG grp;

matrix generate_reference()
{
    matrix ret(9,9,0);
    ret(0,0) = 1.95113485875085;
    ret(0,1) = 0.154783605144974;
    ret(0,2) = 0.0680505878700152;
    ret(0,3) = -0.00794936707126234;
    ret(0,4) = 0.00488252116847035;
    ret(0,5) = -0.016739556907256;
 
    ret(1,1) = 1.37463786228073;
    ret(1,2) = -0.443506551161105;
    ret(1,3) = 0.00434327572533584;
    ret(1,4) = -0.0241623863991684;
    ret(1,5) = -0.0164634137733354;
 
    ret(2,2) = 0.669929538042142;
    ret(2,3) = -0.0178249290507764;
    ret(2,4) = 0.0439660621660147;
    ret(2,5) = 0.0226386579860122;
 
    ret(3,3) = 0.00161865630414999;
    ret(3,4) = -0.000945688539991152;
    ret(3,5) =  -0.000530294431229498;
 
    ret(4,4) = 0.00436245198719005;
    ret(4,5) = 0.0012518352114574;
 
    ret(5,5) = 0.00228388757609715;
 
 
    ret(6,6) = 1.96077638300068;
    ret(6,7) = -0.0960706079156964;
    ret(6,8) = 0.0864761087357964;
 
    ret(7,7) = 1.1910662837907;
    ret(7,8) = 0.305156412769539;
 
    ret(8,8) = 0.844190078267492;

    return ret;
}

matrix compute_diff(MPS<matrix, grp> const & mps, matrix const & ref, std::vector<int> const & site_irreps, std::vector<int> const & config)
{
    size_t L = mps.size();

    matrix rdm(L,L), diff(L,L);
    for (int i=0; i < L-1; ++i)
        for (int j=i+1; j < L; ++j)
        {
            MPO<matrix, grp> mpo = SU2::make_op<matrix, grp>(i, j, site_irreps);
            rdm(i,j) = SU2::expval(mps, mpo, i, j, config);
            diff(i,j) = std::abs(rdm(i,j) - ref(i,j));
        }

    //cout.precision(5);
    //for (int i=0; i < L-1; ++i)
    //{
    //    for (int j=0; j < i+1; ++j) cout << std::setw(10) << " ";
    //    for (int j=i+1; j < L; ++j)
    //    {
    //        cout << std::setw(10) << rdm(i,j);
    //    }
    //    cout << endl;
    //}
    for (int i=0; i < L-1; ++i)
    {
        for (int j=0; j < i+1; ++j) cout << std::setw(10) << " ";
        for (int j=i+1; j < L; ++j)
        {
            cout << std::setw(10) << ref(i,j) / rdm(i,j);
        }
        cout << endl;
    }

    return diff;
}


int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
            return 1;
        }
        MPS<matrix, grp> mps;
        load(argv[1], mps);
        
        size_t L = mps.size();

        std::vector<int> site_irreps;
        for (int i=0; i < L; ++i)
            site_irreps.push_back(mps[i].site_dim()[1].first[2]);
        //std::copy(site_irreps.begin(), site_irreps.end(), std::ostream_iterator<int>(cout, ""));        

        matrix ref = generate_reference();
        std::vector<int> config(6,0);

        for (int i=0; i < L; ++i) {
            MPO<matrix, grp> mpo = SU2::make_count<matrix, grp>(i, site_irreps);
            double n = SU2::expval(mps, mpo, i, i, config);
            maquis::cout <<  n << std::endl;
        }


        //for (int v0=-2; v0 < 3; ++v0) {
        //config[0] = v0;

        //for (int v1=-2; v1 < 3; ++v1) {
        //config[1] = v1;

        //for (int v2=-2; v2 < 3; ++v2) {
        //config[2] = v2;

        //for (int v3=-2; v3 < 3; ++v3) {
        //config[3] = v3;

        //for (int v4=-2; v4 < 3; ++v4) {
        //config[4] = v4;

        //for (int v5=-2; v5 < 3; ++v5) {
        //config[5] = v5;

          matrix diff = compute_diff(mps, ref, site_irreps, config);
          double ss = norm_square(diff);
        //  std::cout << ss << " " << v0 << v1 << v2 << v3 << v4 << v5 << std::endl;
        //}
        //}
        //}
        //}
        //}
        //}
        
    } catch (std::exception& e) {
        std::cerr << "Error:" << std::endl << e.what() << std::endl;
        return 1;
    }
}
