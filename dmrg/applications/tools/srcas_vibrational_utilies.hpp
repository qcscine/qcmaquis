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

#ifndef SAMPLING_VIB_HPP
#define SAMPLING_VIB_HPP

#include <iostream>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <string.h>


void quicksort(string dets[], double b[], int left, int right) {
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

// Calculates CI coefficient for n-mode Determinant via MPS-MPS overlap evaluation
template <typename Determinant, typename Matrix, typename SymmGroup>
double calculate_coefficient(MPS<Matrix, SymmGroup> mps, Determinant det, std::vector<Index<SymmGroup>>& phys_dims, std::vector<int>& site_type, NU1::charge& right_end){
    auto state = HelperClassBasisVectorConverter<NU1>::GenerateIndexFromString(det, phys_dims, site_type, mps.length());
    MPS<Matrix, SymmGroup> mps2 = state_mps<Matrix>(state, phys_dims, site_type, right_end);
    double ovp = overlap(mps, mps2);
    return ovp;
}

// +------------------+
//  SAMPLING STRUCTURE
// +------------------+

class SRCAS {
    // Random number generator, work together with boost library
    // Distribution is a variable modelling a random distribution over the 0-1 range
    SRCAS() : distribution(0.,1.), random_number(generator,distribution)
    {
        unsigned int seed = 123456;
        generator.seed(seed+time(NULL));
    }

    // Definition of the variables used for the random generation
    boost::mt19937 generator;
    boost::uniform_real<> distribution;
    boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > random_number;

    // -- ARGUMENTS of perform_srcas --
    // 1)  mps           --> reference MPS
    // 2)  modals        --> modal basis size for each mode, provided in input in determinant format
    // 3)  starting_det  --> starting determinant provided in input
    // 4)  nsample       --> number of samples per macroiteration
    // 5)  nitermax      --> maximum number of macroiterations of the algorithm
    // 6) CI_threshold   --> threshold to use to store the determinant
    // 7) COM_threshold  --> threshold to assess the convergence of the SRCAS algorithm

    /*
    add_option("init_basis_state", "local indexes for basis state init", value(""));

         // Vibrational SRCAS
        add_option("srcas_targetCompleteness", "Desired completness for SRCAS to terminate sampling", value(0.9));
        add_option("srcas_maxNumIterations", "Maximum number of macroiterations until SRCAS sampling is terminated", value(10));
        add_option("srcas_numMicroiterations", "Number of microiterations in each SRCAS macroiteration", value(10000));
    */

    template <typename Determinant, typename Matrix, typename SymmGroup>
    void perform_srcas(MPS<Matrix, SymmGroup> mps, std::vector<int> modals, Determinant starting_det, int nsamples, int nitermax, double CI_threshold, double COM_threshold)
    { 
        // Type definition
        typedef typename MPS<Matrix,SymmGroup>::scalar_type scalar_type ;
        // Determinants initialization
        std::size_t nmodes = starting_det.size();
        Determinant det_queen, det_tmp;
        det_queen.resize(nmodes);
        det_tmp.resize(nmodes);
        det_queen = starting_det;
        // This is used for determinants list
        typedef std::map< Determinant , double> Hash_Map_with_value;
        Hash_Map_with_value hash;
        typename Hash_Map_with_value::iterator iter;

        // Stuff that is needed to construct mps for each det
        std::vector<int> site_type; // this vector contains the mode to which each modal belongs to, e.g. 0,0, ... ,0,1,1, ... for an nmode MPS
        std::vector<int> mode_start; // this vector contains the indices of the first modal of each mode on the lattice
        for (std::size_t idx1 = 0; idx1 < modals.size(); idx1++) {
            if (idx1 == 0) mode_start.push_back(0);
            else mode_start.push_back(mode_start[idx1-1]+modals[idx1-1]);
            for (std::size_t idx2 = 0; idx2 < modals[idx1]; idx2++) {
                site_type.push_back(static_cast<int>(idx1));
            }
        }
        std::vector<Index<SymmGroup>> phys_dims;
        for (std::vector<int>::iterator it = mode_start.begin(); it != mode_start.end(); ++it){
            phys_dims.push_back(mps[*it].site_dim());
        }
        NU1::charge right_end = NU1::charge(std::vector<int>(phys_dims.size(), 1));
    
        // Initialize variables that are used during the sampling
        scalar_type ci, ci0, ci_ratio, ci_tmp;
        scalar_type sum_ci2 = 0.0, uncompleteness = 0.0;
        int nMacroIter = 0, nsampled = 0, naccepted_queen = 0;
        int nmodes_excited, mode_2excite;
        float x;
        // +-----------+
        //   MAIN LOOP
        // +-----------+
        do {
            // For every macroiteration generate N determinants
            for ( int isample = 0; isample < nsamples ; isample++ ) {
                // Loop over the determinants to be sampled
                nmodes_excited = int(ceil(nmodes*random_number())); // randomly pick how many modes will be (de-)excited
                det_tmp = det_queen;
                for (std::size_t idx = 0; idx < nmodes_excited; idx++){
                    // Compute mode excitation
                    mode_2excite = int(floor(nmodes*random_number())); // pick a random mode
                    det_tmp[mode_2excite] = int(floor(modals[mode_2excite]*random_number())); // pick a random modal for that mode --> could be changed to be closer to queen with only +-1
                } ;
                // Updates the data if the determinant has not been visited yet.
                iter = hash.find(det_tmp) ;
                if(iter == hash.end()) {
                    ci = calculate_coefficient<Determinant,Matrix,SymmGroup>(mps, det_tmp, phys_dims, site_type, right_end);
                    // The data are stored based on the CI_threshold parameter
                    if(std::fabs(ci) >= CI_threshold) {
                        hash[det_tmp] = ci ;
                        nsampled++;
                    }
                } else {
                    ci = iter->second;
                }
                // Determinant update (regardless of being in the hash table or not to avoid getting stuck in the Markov chain)
                ci_ratio = pow(ci,2.0)/pow(ci0,2) ;
                x = random_number() ;
                if( ci_ratio > x ) {
                    det_queen = det_tmp ;
                    ci0 = ci ;
                    naccepted_queen++ ;
                }
            }
            sum_ci2 = 0.0 ;
            for (iter=hash.begin(); iter!=hash.end(); iter++) {
                ci_tmp  = iter->second;
                sum_ci2 += pow(ci_tmp,2.0);
            }
            nMacroIter++ ;
            // Prints results
            maquis::cout << "----------------------------------------------------------------" << std::endl ;
            maquis::cout << "Macroiteration number:                       " << nMacroIter << std::endl;
            maquis::cout << "Determinants sampled above the CI threshold: " << nsampled << std::endl;
            maquis::cout << "Determinants accepted as queens:             " << naccepted_queen << std::endl;
            maquis::cout << "Current completeness (\\sum(ci^2)):           " << sum_ci2 << std::endl;
            uncompleteness = 1.0-sum_ci2;
        } while( uncompleteness > COM_threshold && nMacroIter< nitermax ) ;

        

        // +---------------+
        //   FINAL PRINTING
        // +---------------+
        double    CIs_show[hash.size()]; // CI value
        string    dets_show[hash.size()]; // dets represent
        int i = 0;
        int det_length = hash.begin()->first.size();
        maquis::cout << std::endl << "-------------DETERMINANTS ABOVE CI THRESHOLD--------------------" << std::endl << std::endl;
        for(iter = hash.begin(); iter!=hash.end(); iter++) {
            // Local initialization that is later used for sorting
            string ctmp;
            CIs_show[i] = iter->second;
            for(int p = 0; p < det_length; p++)
                ctmp = ctmp + boost::lexical_cast<string>(iter->first[p]);
            dets_show[i] = ctmp;
            i++;
        }
        // Final sorting
        quicksort(dets_show, CIs_show,0, hash.size()-1);
        // Output the entire ordered list
        maquis::cout << std::fixed << std::setprecision(10);
        for(int i = 0; i < hash.size() ; i++){
            maquis::cout << " Determinant " << dets_show[hash.size()-i-1] << " with ";
            if (CIs_show[hash.size()-i-1]>0) maquis::cout << " ";
            maquis::cout << CIs_show[hash.size()-i-1] << " is number " << i+1 << std::endl;
        }
        
    }
};



#endif


