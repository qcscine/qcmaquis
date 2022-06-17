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

#include "srcas_utilities.h"
#include "maquis_dmrg.h"


#include <iostream>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <string>

// Random number generator, work together with boost library
// Distribution is a variable modelling a random distribution over the 0-1 range
SRCAS::SRCAS(DmrgParameters& parameters) : distribution_(0.,1.), randomNumber_(generator_,distribution_), parms_(parameters)
{
    generator_.seed(parms_["seed"]);
    // Get the number of modes
    if(parms_["MODEL"] == "nmode")
        numModes_ = parms_["nmode_num_modes"];
    else if (parms_["MODEL"] == "watson")
        numModes_ = parms_["L"];
     // If user set a starting det, use this, otherwise use "0,0,0, ... ,0"
    if (parms_.is_set("init_basis_state"))
        startingDet_ = parms_["init_basis_state"].str();
    else {
        startingDet_ = "0";
        for (int i=1; i<numModes_; i++) startingDet_ += ",0";
    }
    detQueen_.resize(numModes_);
    detTmp_.resize(numModes_);
    for (int i=0; i<detQueen_.size(); i++) {
        std::cout << startingDet_[2*i] << " ";
        detQueen_[i]=startingDet_[2*i];
        std::cout << detQueen_[i] << std::endl;
    }
}

void SRCAS::printSRCASSettings() {
    std::cout << "--- SRCAS SETTINGS ---" << std::endl;
    std::cout << "MPS taken from:                      " << parms_["chkpfile"].str() << std::endl;
    std::cout << "Determinant space is:                ";
    if(parms_["MODEL"] == "nmode") std::cout << parms_["nmode_num_basis"].str() << std::endl;
    else if(parms_["MODEL"] == "watson")  {
        std::cout << parms_["Nmax"];
        for (int i=1; i<numModes_; i++) std::cout << "," << parms_["Nmax"];
        std::cout << std::endl;
    }
    std::cout << "Starting determinant is:             " << startingDet_ << std::endl;
    std::cout << "SRCAS target completness is:         " << std::setprecision(2) << std::fixed << parms_["srcas_targetCompleteness"] << std::endl;
    std::cout << "Maximum number of iterations is:     " << parms_["srcas_maxNumIterations"] << std::endl;
    std::cout << "Number of samples per iteration is:  " << parms_["srcas_numSamples"] << std::endl;
    std::cout << "Random number seed is:               " << parms_["seed"] << std::endl;
}

void SRCAS::quicksort(std::string dets[], double b[], int left, int right) {
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
            std::string ctmp ;
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



void SRCAS::run() {
    std::cout << std::endl << "-- Starting SRCAS --" << std::endl;
    // Creates the interface object
    maquis::DMRGInterface<double> interface(parms_);

    auto overlap = interface.getCICoefficient(startingDet_);
    std::cout << overlap << std::endl;

    // Initialize variables that are used during the sampling
    //scalar_type ci, ci0, ci_ratio, ci_tmp;
    //scalar_type sum_ci2 = 0.0, uncompleteness = 0.0;
    double completeness;
    int nMacroIter = 0, nSampled = 0, nAcceptedQueen = 0;
    

    // +-----------+
    //   MAIN LOOP
    // +-----------+
    do {
        // For every macroiteration generate N determinants
        for ( int isample = 0; isample < parms_["srcas_numSamples"] ; isample++ ) {
            // Loop over the determinants to be sampled
            int numModesToExcite = int(ceil(numModes_*randomNumber_())); // randomly pick how many modes will be (de-)excited
            std::cout << numModesToExcite << " ";
            
            detTmp_= detQueen_;
            for (std::size_t idx = 0; idx < numModesToExcite; idx++){
                // Compute mode excitation
                int modeToExcite = int(floor(numModes_*randomNumber_())); // pick a random mode
                std::cout << modeToExcite << " ";
                //detTmp_[modeToExcite] = int(floor(modals[modeToExcite]*randomNumber_())); // pick a random modal for that mode --> could be changed to be closer to queen with only +-1
            } ;
            /*
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
            */
        }
        /*
        sum_ci2 = 0.0 ;
        for (iter=hash.begin(); iter!=hash.end(); iter++) {
            ci_tmp  = iter->second;
            sum_ci2 += pow(ci_tmp,2.0);
        }
        */
        nMacroIter++ ;
        // Prints results
        maquis::cout << "----------------------------------------------------------------" << std::endl ;
        maquis::cout << "Macroiteration number:                       " << nMacroIter << std::endl;
        maquis::cout << "Determinants sampled above the CI threshold: " << nSampled << std::endl;
        maquis::cout << "Determinants accepted as queens:             " << nAcceptedQueen << std::endl;
        //maquis::cout << "Current completeness (\\sum(ci^2)):          " << completeness << std::endl;
        
    } while((completeness < parms_["srcas_targetCompleteness"]) && (nMacroIter < parms_["srcas_maxNumIterations"]));

/*
        // Opens the determinant file and loops over it
        std::string nameOfDetFile = opt.parms["determinant_file"];
        double threshold = opt.parms["determinant_threshold"];
        std::ifstream is(nameOfDetFile);
        std::string str;
        while (getline(is, str)) {
            
            if (std::abs(overlap) > threshold)
                std::cout << "CI coefficient of " << str << " : " << overlap << std::endl;
        }
    
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
*/
}    

        // +---------------+
        //   FINAL PRINTING
        // +---------------+
/*
void SRCAS::printResults() {
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
*/
