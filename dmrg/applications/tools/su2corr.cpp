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

#include <boost/lambda/lambda.hpp>

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

std::vector<double> generate_2rdm_ref()
{
    // reference for <cdag_i cdag_i+1 c_i+2 c_i+3>
    std::vector<double> ret;
    //ret.push_back(0.0057801190338182);      //0123
    //ret.push_back(0.000504560118115581);    //1234
    //ret.push_back(3.94702812535319e-05);    //2345
    //ret.push_back(0);                       //3456
    //ret.push_back(0.000110219364671496);    //4567
    //ret.push_back(0);                       //5678

    // reference from Maquis
    ret.push_back(-8.623925e-03);  //0123
    ret.push_back(2.521896e-04);   //1234
    ret.push_back(1.085327e-06);   //2345
    ret.push_back(0);              //3456
    ret.push_back(-5.364513e-04);  //4567
    ret.push_back(0);              //5678

    return ret;
}

template <class Matrix>
void print_triang(Matrix const & mat)
{
    size_t L = num_cols(mat);
    for (int i=0; i < L-1; ++i)
    {
        for (int j=0; j < i+1; ++j) cout << std::setw(10) << " ";
        for (int j=i+1; j < L; ++j)
        {
            cout << std::setw(10) << mat(i,j);
        }
        cout << endl;
    }
}

template <class Matrix, class SymmGroup>
matrix compute_ratio(MPS<Matrix, SymmGroup> const & mps, matrix const & ref, std::vector<int> const & site_irreps, std::vector<int> const & config)
{
    size_t L = mps.size();

    Matrix rdm(L,L), ratio(L,L);
    for (int i=0; i < L-1; ++i)
        for (int j=i+1; j < L; ++j)
        {
            MPO<Matrix, SymmGroup> mpo = SU2::make_1rdm_term<Matrix, SymmGroup>(i, j, site_irreps);
            rdm(i,j) = SU2::expval(mps, mpo, i, j, config);
            ratio(i,j) = ref(i,j) / rdm(i,j);
        }

    return ratio;
}

template <class Matrix, class SymmGroup>
std::vector<typename Matrix::value_type>
compute_off_diag_ratio(MPS<Matrix, SymmGroup> const & mps, int which_diag, matrix const & ref,
                              std::vector<int> const & site_irreps, std::vector<int> const & config)
{
    size_t L = mps.size() - which_diag;
    std::vector<typename Matrix::value_type> ret(L);
    for (int i=0; i < L; ++i)
    {
        MPO<Matrix, SymmGroup> mpo = SU2::make_1rdm_term<Matrix, SymmGroup>(i, i+which_diag, site_irreps);
        ret[i] = SU2::expval(mps, mpo, i, i+which_diag, config);
        ret[i] = ref(i,i+which_diag) / ret[i];
    }
    return ret;
}

int main(int argc, char ** argv)
{
    try {
        if (argc != 2) {
            std::cout << "Usage: " << argv[0] << " <mps.h5>" << std::endl;
            return 1;
        }
        //cout.precision(5);

        MPS<matrix, grp> mps;
        load(argv[1], mps);
        size_t L = mps.size();

        std::vector<int> site_irreps;
        for (int i=0; i < L; ++i)
            site_irreps.push_back(mps[i].site_dim()[1].first[2]);

        matrix ref = generate_reference();
        std::vector<int> config(20,0);

        //for (int i=0; i < L; ++i) {
        //    MPO<matrix, grp> mpo = SU2::make_count<matrix, grp>(i, site_irreps);
        //    double n = SU2::expval(mps, mpo, i, i, config);
        //    maquis::cout <<  n << std::endl;
        //}


        matrix ratio = compute_ratio(mps, ref, site_irreps, config);
        print_triang(ratio);

        //int i=0, j=3;
        //MPO<matrix, grp> mpo = SU2::make_custom<matrix, grp>(i, j, site_irreps);
        //double eval = SU2::expval(mps, mpo, i, j, config);
        //maquis::cout << eval << "  ref: " << ref(i,j) + ref(i,j+1) << std::endl;

        //std::vector<double> od = compute_off_diag_ratio(mps, 2, ref, site_irreps, config);
        //std::copy(od.begin(), od.end(), std::ostream_iterator<double>(cout, "  "));
        //std::transform(od.begin(), od.end(), std::ostream_iterator<double>(cout, "  "), boost::lambda::_1 / od[0]);
        //cout << endl;

        //for (int v1=0; v1 < 2; ++v1) {
        //for (int v2=0; v2 < 2; ++v2) {

        //for (int v3=0; v3 < 8; ++v3) {
        //for (int v4=0; v4 < 4; ++v4) {
        //for (int v5=0; v5 < 4; ++v5) {

        //std::vector<double> ref2 = generate_2rdm_ref();
        //std::vector<double> result;
        //for(int i=0; i < L-3; ++i)
        //{
        //    //MPO<matrix, grp> four = SU2::make_2rdm_term_custom<matrix, grp>(i,i+1,i+2,i+3, 1,1,v3,v4,v5, site_irreps);
        //    MPO<matrix, grp> four = SU2::make_2rdm_term_custom<matrix, grp>(i,i+1,i+2,i+3, 1,1,0,0,0, site_irreps);
        //    double twodm0123 = SU2::expval(mps, four, i+10,0, config);
        //    result.push_back(twodm0123 / ref2[i]);
        //}
        ////std::copy(result.begin(), result.end(), std::ostream_iterator<double>(cout, "  "));
        //std::transform(result.begin(), result.end(), std::ostream_iterator<double>(cout, "  "), boost::lambda::_1 / result[0]);
        //std::cout << std::endl;

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
