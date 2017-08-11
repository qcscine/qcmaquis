/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Yingjin Ma <yingjin.ma@phys.chem.ethz.ch>
 *               2017 by Alberto Baiardi <alberto.baiardi@sns.it>
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

#ifndef SAMPLING_HPP
#define SAMPLING_HPP

#include <iostream>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <string.h>

#include "ci_encode.hpp"

using std::cerr;
using std::cout;
using std::endl;
using std::string;

// refer from wiki
void quick_sort (double data[], int left, int right) {
    // Initialization
	int p = (left + right) / 2;
	int pivot = data[p];
	int i = left, j = right ;
	for ( ; i < j;) {
		while (! (i>= p || std::abs(pivot) < std::abs(data[i])))
			++i;
		if (i < p) {
			data[p] = data[i];
			p = i;
		}
		while (! (j <= p || std::abs(data[j]) < std::abs(pivot)))
			--j;
		if (j > p) {
			data[p] = data[j];
			p = j;
		}
	}
	data[p] = pivot;
	if (p - left > 1)
		quick_sort(data, left, p - 1);
	if (right - p > 1)
		quick_sort(data, p + 1, right);
}

void quicksort(string dets[], double b[], int left, int right) {
    //
    double pivot = std::abs(b[(left+right)/2]);
    int l = left;
    int r = right;
    while(l <= r) {
        while (std::abs(b[l]) < pivot)
            l++ ;
        while (std::abs(b[r]) > pivot)
            r-- ;
        if (l <= r) {
            // Variable definition
            double tmp ;
            tmp  = b[l] ;
            b[l] = b[r] ;
            b[r] = tmp ;
            string ctmp ;
            ctmp = dets[l] ;
            dets[l] = dets[r] ;
            dets[r] = ctmp ;
            l++ ;
            r-- ;
        }
    };
    // Calls the routine defined above
    if (left < r)  quicksort(dets, b, left, r);
    if (l < right) quicksort(dets, b, l, right);
}

// +------------------+
//  SAMPLING STRUCTURE
// +------------------+

struct Sampling {
    // Random number generator, work together with boost library
    // Distribution is a variable modelling a random distribution over the 0-1 range
    Sampling(int id) : distribution(0.,1.), random_number(generator,distribution)
    {
        unsigned int seed = 123456;
        generator.seed(seed/(id+1)+time(NULL));
    }
    // No integer provided in input
    Sampling() : distribution(0.,1.), random_number(generator,distribution)
    {
        unsigned int seed = 123456;
        generator.seed(seed+time(NULL));
    }
    // template generate_dets
    template <typename Determinants, typename Matrix, typename SymmGroup, typename Hash_value, typename Hash_index>
    void generate_dets(Determinants dets, Determinants dets_mclr, MPS<Matrix, SymmGroup> mps,
                       Hash_value hash, Hash_index hash_index, std::vector< Index<SymmGroup> > site_dims,
                       int norb, int nsample, double CI_threshold, double COM_threshold)
    {
        // Header printing
        maquis::cout << "    CI-threshold  : " <<  CI_threshold << std::endl;
        maquis::cout << "   COM-threshold  : " << COM_threshold << std::endl;
        // Type definition
        typename Hash_value::iterator                       iter ;
        typename Hash_index::iterator                       iter_idx ;
        typedef typename Determinants::value_type           Determinant ;
        typedef typename MPS<Matrix,SymmGroup>::scalar_type scalar_type ;
        // Set up the initial parameters, fix the symmetry group
        scalar_type ci , ci0 , ci_ratio , ci_tmp ;
        scalar_type sum_ci2 = 0.0 ;
        scalar_type completeness = 0.0 ;
        std::size_t det_length = dets[0].size() ;
        std::size_t n_total_dim ;
	    std::size_t number_of_dets ;
        // This part will be used as the determinants reservoir
        Determinant det;
        for (std::size_t c = 0; c < dets.size(); ++c)
            hash_index[dets[c]] = c ;
        // determinant spawnning -- preparing part
        Determinant det_queen , det_bee , det_tmp ; // determinant "queen bee tmp"
        int iaccept = 0 ;
        int iaccept_queen = 0 ;
        // Loop over the determinants to be followed
        for (std::size_t c = 0; c < dets_mclr.size(); ++c) {
            det = dets_mclr[c];
            //ci0 = extract_coefficient(mps, det);
            hash[det] = ci0 ;
        }
        // Computes the total number of electron/holes for alpha or beta orbitals
        n_total_dim = 0 ;
        for (std::size_t i = 0; i < det_length; i++)
            n_total_dim += site_dims[i][0].second ;
        // Get the number of excited electrons. This number is generated randomly
        int nmodes_excited ;
        det_queen = dets[0];
        int nMAX = 0 ;
        // +-----------+
        //   MAIN LOOP
        // +-----------+
        do {
            nMAX++ ;
            for ( int isample = 0; isample < nsample; isample++ ) {
                // Loop over the determinants to be sampled
                // Updates the "queen" determinant
                nmodes_excited = int(floor(n_total_dim*random_number()) + 1);
                det_tmp = det_queen ;
                int iele_excited = 0 ;
                for (std::size_t idx = 0; idx < nmodes_excited; idx++){
                    // Excites an alpha electron with probability 1/2
                    iele_excited++ ;
                    // Decides the electron to be excited
                    int mode_2excite = int(floor(random_number()*n_total_dim )) ;
                    det_tmp[mode_2excite] += 1 ;
                } ;
                // Corrects the occupation number by the modulus
                for (std::size_t idx = 0; idx < det_length; idx++) {
                    det_tmp[idx] %= site_dims[idx][0].second ;
                    std::cout << det_tmp[idx] << std::endl ;
                }
                // Updates the data
                iter = hash.find(det_tmp) ;
                std::cout << "Pippo" << std::endl ;
                if(iter == hash.end()) {
                    ci = extract_coefficient<Matrix,SymmGroup>(mps, det_tmp);
                    std::cout << ci << std::endl ;
                    std::cout << "Pippo" << std::endl ;
                    if(std::abs(ci) >= CI_threshold) {
                        std::cout << "Pippo" << std::endl ;
                        hash[det_bee] = ci ;
                        iaccept++;
                    }
                } else {
                    ci = iter->second;
                }
                sum_ci2  = sum_ci2 + pow(ci,2.0);
                ci_ratio = pow(ci,2.0)/pow(ci0,2);
                // Whether use this bee-det as the new queen-det
                if( ci_ratio > random_number()) {
                    det_queen = det_tmp ;
                    ci0 = ci ;
                    iaccept_queen++ ;
                }
            }
            sum_ci2=0.0;
            for( iter=hash.begin(); iter!=hash.end(); iter++) {
                ci_tmp  = iter->second;
                sum_ci2 = sum_ci2+pow(ci_tmp,2.0);
            }
            // Prints results
            maquis::cout << "Determinant-naccept" << iaccept << std::endl;
            maquis::cout << "Determinant-naccept-queen" << iaccept_queen << std::endl;
            completeness=1.0-sum_ci2;
            maquis::cout << " Current completeness (1-\\sum(ci^2)) : " << completeness << std::endl;
        } while(completeness>COM_threshold && nMAX< 100) ;
        // +---------------+
        //   FINAL PRINTING
        // +---------------+
        double    CIs_show[hash.size()]  ; //value
        uint32_t  CIs_index[hash.size()] ; //index in the reservoir
        string    dets_show[hash.size()] ; //dets represent
        uint32_t  i = 0 ;
        for(iter=hash.begin(); iter!=hash.end(); iter++) {
            // Local initialization
            string ctmp;
            CIs_show[i] = iter->second;
            for(int p = 0; p < det_length; p++)
                ctmp = ctmp + boost::lexical_cast<string>(iter->first[p]) ;
            dets_show[i] = ctmp;
            iter_idx = hash_index.find(iter->first);
            if(iter_idx == hash_index.end()) {
                CIs_index[i] = 0 ;
            } else {
                CIs_index[i] = iter_idx->second + 1 ;
            }
            i++;
        }
        std::ofstream fout;
        fout.open("det_coeff.tmp");
        fout << i << std::endl;
        for( uint32_t i=0 ; i<hash.size() ; i++)
            if(CIs_index[hash.size()-i-1]!=0)
                fout << CIs_index[hash.size()-i-1] << "  " << CIs_show[hash.size()-i-1] << std::endl;
        // Final sorting
        quicksort(dets_show, CIs_show,0, hash.size()-1);
        for(int i=0; i<hash.size() ; i++)
            maquis::cout << " The sorted determinants: " << dets_show[hash.size()-i-1] << " " << CIs_show[hash.size()-i-1] << std::endl;
    }
    // Definition of the variables used for the random generation
    boost::mt19937 generator;
    boost::uniform_real<> distribution;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<> > random_number;
};



#endif


