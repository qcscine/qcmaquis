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

#ifndef SAMPLING_VIB_H
#define SAMPLING_VIB_H

#include "dmrg/utils/DmrgParameters.h"

#include <string>
#include <boost/random.hpp>

class SRCAS {
    public:
        SRCAS(DmrgParameters& parameters);      
        void run();
        void printSRCASSettings();
        void printResults();
    private:
        void quicksort(std::string dets[], double b[], int left, int right);

        boost::mt19937 generator_;
        boost::uniform_real<> uniformDist_;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > uniformRandomNumber_;

        DmrgParameters& parms_;
        std::string startingDet_, maxDetStr_, detTmpStr_;
        std::vector<int> detQueen_, detTmp_, detSpace_;
        int numModes_;
        double completeness_;
        const double samplingFraction_ = 1.0/3.0; // Change this value to alter the speed of sampling across the Hilbert space
        std::map<std::vector<int>, double> hashTable_;
        std::map<std::vector<int>, double>::iterator iter_;
};

#endif