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
#include "dmrg/utils/DmrgParameters.h"

#include <iostream>
#include <iomanip>
#include <boost/random.hpp>
#include <boost/lexical_cast.hpp>
#include <math.h>
#include <string>

template <typename ScalarType> // real or complex
SRCAS<ScalarType>::SRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface) :
    interface_(interface), uniformDist_(0.,1.), uniformRandomNumber_(generator_,uniformDist_), parms_(parameters)
{
    generator_.seed(parms_["seed"]);
    // Get the number of modes and the maximum occupation of each one
    if(parms_["MODEL"] == "nmode") {
        numModes_ = parms_["nmode_num_modes"];
        maxDetStr_ = parms_["nmode_num_basis"].str();
        detSpace_ = parms_["nmode_num_basis"].as<std::vector<int> >();
    } else if (parms_["MODEL"] == "watson") {
        numModes_ = parms_["L"];
        maxDetStr_ = parms_["Nmax"].str();
        detSpace_ = parms_["Nmax"].as<std::vector<int> >();
        if (detSpace_.size()!=numModes_ && detSpace_.size()!=1){
            throw std::runtime_error("The Nmax parameter must be either a single integer, or a vector of lenght L");
        }
        if (detSpace_.size()!=numModes_) {
            for (int i=1; i<numModes_; i++) {
                maxDetStr_ += ",";
                maxDetStr_ += parms_["Nmax"].str();
            }
            std::vector<int> tmpVec(numModes_, std::stoi(parms_["Nmax"].str()));
            detSpace_ = std::move(tmpVec);
        }
    }

    // If user set a starting det, use this, otherwise use "0,0,0, ... ,0"
    if (parms_.is_set("init_basis_state")) {
        startingDet_ = parms_["init_basis_state"].str();
        detQueen_ = parms_["init_basis_state"].as<std::vector<int> >();
    } else {
        startingDet_ = "0";
        for (int i=1; i<numModes_; i++) startingDet_ += ",0";
        std::vector<int> tmpVec(numModes_, 0);
        detQueen_ = std::move(tmpVec);
    }
    detTmp_ = detQueen_;
}

template <typename ScalarType> // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::printSRCASSettings() {
    maquis::cout << std::endl << "----- SRCAS SETTINGS -----" << std::endl;
    maquis::cout << "MPS taken from:                            " << parms_["chkpfile"].str() << std::endl;
    maquis::cout << "Determinant space is:                      " << maxDetStr_ << std::endl;
    maquis::cout << "Starting determinant is:                   " << startingDet_ << std::endl;
    maquis::cout << "CI coeff (overlap) threshold is:           " << parms_["srcas_overlapThreshold"] << std::endl;
    maquis::cout << "SRCAS target completeness is:              " << parms_["srcas_targetCompleteness"] << std::endl;
    maquis::cout << "Maximum number of iterations is:           " << parms_["srcas_maxNumIterations"] << std::endl;
    maquis::cout << "Number of samples per iteration is:        " << parms_["srcas_numSamples"] << std::endl;
    maquis::cout << "Random number seed is:                     " << parms_["seed"] << std::endl;
    maquis::cout << "Fraction of proposed updates to accept is: " << parms_["srcas_samplingFraction"] << std::endl;
}

template <typename ScalarType> // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::quicksort(std::string dets[], ScalarType b[], int left, int right) {
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
            ScalarType tmp ;
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

template <typename ScalarType> // real or complex, nmode or canonical (watson)
std::vector<int> SRCAS<ScalarType>::generateNewDet() {
    // Start from queen
    detTmp_= detQueen_;
    if(parms_["MODEL"] == "nmode" || parms_["MODEL"] == "watson") {
        // Loop over the modes            
        for (int i=0; i<detTmp_.size(); i++) {
            boost::poisson_distribution<> poissonDist(detTmp_[i]+0.5); // poisson distribution centered on the current modal
            boost::variate_generator<boost::mt19937&, boost::poisson_distribution<>> poissonRandomNumber(generator_,poissonDist);
            do {
                if (uniformRandomNumber_() < parms_["srcas_samplingFraction"]) // Only accept a fraction of the proposed updates to stay closer to reference det
                    detTmp_[i] = poissonRandomNumber();
            } while (!(detTmp_[i] < detSpace_[i])); // Only accept valid occupations
        }
    } else {
        maquis::cout << "SRCAS determinant generation NYI for non-vibrational calculations! Abort!" << std::endl;
        exit(1);
    }
    return detTmp_;
}


template <typename ScalarType> // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::run() {
    maquis::cout << std::endl << "----- Starting SRCAS -----" << std::endl << std::endl;

    // Starting det should always be added to the list
    ScalarType overlap = interface_->getCICoefficient(startingDet_);
    hashTable_[detQueen_] = overlap;

    // Initialize variables that are used during the sampling
    double x, ci_ratio, sum_ci2 = 0.0;
    ScalarType ci_tmp, ci0 = overlap;
     
    int nMacroIter = 0, nSampled = 1, nAcceptedQueen = 0;    

    // +-----------+
    //   MAIN LOOP
    // +-----------+
    do {
        // For every macroiteration generate N determinants
        for (int isample = 0; isample < parms_["srcas_numSamples"]; isample++) {
            // Get new determinant
            detTmp_ = generateNewDet();

            // Updates the data if the determinant has not been visited yet.
            iter_ = hashTable_.find(detTmp_) ;
            if(iter_ == hashTable_.end()) {
                detTmpStr_ = std::to_string(detTmp_[0]);
                for (int i=1; i<detTmp_.size(); i++) {
                    detTmpStr_ += ",";
                    detTmpStr_ += std::to_string(detTmp_[i]);
                }
                overlap = interface_->getCICoefficient(detTmpStr_);
                // The data are stored based on the CI_threshold parameter
                if(std::abs(overlap) >= parms_["srcas_overlapThreshold"]) {
                    hashTable_[detTmp_] = overlap;
                    nSampled++;
                }
            } else {
                overlap = iter_->second;
            }
            // Determinant update (regardless of being in the hash table or not to avoid getting stuck in the Markov chain)
            // Selection criterion based on CI coeff^2 in analogy to the completeness measure
            ci_ratio = pow(std::abs(overlap),2.0)/pow(std::abs(ci0),2);
            x = uniformRandomNumber_();
            if (ci_ratio > x) {
                detQueen_ = detTmp_;
                ci0 = overlap;
                nAcceptedQueen++ ;
            }
        }
        sum_ci2 = 0.0 ;
        for (iter_=hashTable_.begin(); iter_!=hashTable_.end(); iter_++) {
            ci_tmp  = iter_->second;
            sum_ci2 += pow(std::abs(ci_tmp),2.0);
        }
        nMacroIter++ ;
        
        // Prints results
        maquis::cout << "----------------------------------------------------------------" << std::endl ;
        maquis::cout << "Macroiteration number:                       " << nMacroIter << std::endl;
        maquis::cout << "Determinants sampled above the CI threshold: " << nSampled << std::endl;
        maquis::cout << "Determinants accepted as queens:             " << nAcceptedQueen << std::endl;
        maquis::cout << "Current completeness (\\sum(ci^2)):           " << sum_ci2 << std::endl;
        
    } while((sum_ci2 < parms_["srcas_targetCompleteness"]) && (nMacroIter < parms_["srcas_maxNumIterations"]));
    // Final completeness
    completeness_ = sum_ci2;
}    

// +---------------+
//   FINAL PRINTING
// +---------------+
template <typename ScalarType> // real or complex, nmode or canonical (watson)
void SRCAS<ScalarType>::printResults() {
    maquis::cout << "----------------------------------------------------------------" << std::endl ;
    maquis::cout << std::endl << "--- Finished SRCAS ---" << std::endl;
    maquis::cout << "Final completeness is:                " << completeness_ << std::endl;
    maquis::cout << "# of stored determinants is:          " << hashTable_.size() << std::endl;

    ScalarType CIs_show[hashTable_.size()]; // CI value
    std::string dets_show[hashTable_.size()]; // dets represent
    int i = 0;
    int det_length = hashTable_.begin()->first.size();
    maquis::cout << std::endl << "-------------DETERMINANTS ABOVE OVERLAP THRESHOLD--------------------" << std::endl << std::endl;
    for(iter_ = hashTable_.begin(); iter_!=hashTable_.end(); iter_++) {
        // Local initialization that is later used for sorting
        std::string ctmp;
        CIs_show[i] = iter_->second;
        for(int p = 0; p < det_length; p++)
            ctmp = ctmp + boost::lexical_cast<std::string>(iter_->first[p]);
        dets_show[i] = ctmp;
        i++;
    }
    // Final sorting
    quicksort(dets_show, CIs_show, 0, hashTable_.size()-1);
    // Output the entire ordered list
    maquis::cout << std::fixed << std::setprecision(10);
    for(int i = 0; i < hashTable_.size() ; i++){
        maquis::cout << " Determinant " << dets_show[hashTable_.size()-i-1] << " with ";
        //if (CIs_show[hashTable_.size()-i-1]>0) maquis::cout << " ";
        maquis::cout << CIs_show[hashTable_.size()-i-1] << " is number " << i+1 << std::endl;
    }        
}

template <typename ScalarType> // real or complex, nmode or canonical (watson)
std::vector<int> SRCAS<ScalarType>::getCurrentQueen() {
    return detQueen_;
}

template <typename ScalarType> // real or complex, nmode or canonical (watson)
std::map<std::vector<int>, ScalarType> SRCAS<ScalarType>::getDetTable() {
    return hashTable_;
}

template <typename ScalarType> // real or complex, nmode or canonical (watson)
double SRCAS<ScalarType>::getCompleteness() {
    return completeness_;
}

// Explicit template instantiation
template class SRCAS<double>;
template class SRCAS<std::complex<double>>;