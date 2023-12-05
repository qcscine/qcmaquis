/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef SAMPLING_VIB_H
#define SAMPLING_VIB_H

#include "dmrg/utils/DmrgParameters.h"
#include "maquis_dmrg.h"

#include <string>
#include <memory>
#include <boost/random.hpp>

template <typename ScalarType> // real or complex
class SRCAS {
    using InterfaceType = maquis::DMRGInterface<ScalarType>;
    public:
        SRCAS(DmrgParameters& parameters, std::shared_ptr<InterfaceType> interface);      
        void run();
        void printSRCASSettings();
        void printResults();

        std::vector<int> getCurrentQueen();
        std::map<std::vector<int>, ScalarType> getDetTable();
        double getCompleteness();

    private:
        std::vector<int> generateNewDet();
        double calculateCompleteness();
        double addToCompleteness (ScalarType coeff, std::vector<int> det);
        void quicksort(std::string dets[], ScalarType b[], int left, int right);
        
        // For electronic case
        int getARandomOccSpinOrb(std::vector<int> det);
        int getARandomUnoccSpinOrb(std::vector<int> det);
        bool symmetriesFulfilled(std::vector<int> det);

        boost::mt19937 generator_;
        boost::uniform_real<> uniformDist_;
        boost::geometric_distribution<double> geomDist_;
        boost::variate_generator<boost::mt19937&, boost::uniform_real<double> > uniformRandomNumber_;
        boost::variate_generator<boost::mt19937&, boost::geometric_distribution<double> > geometricRandomNumber_;

        DmrgParameters& parms_;
        std::shared_ptr<InterfaceType> interface_;

        std::string startingDet_, maxDetStr_, detTmpStr_;
        std::vector<int> detQueen_, detTmp_, detSpace_;
        int numParticles_; // this is either modes or electrons, for the vibrational or the electronic case, respectively
        double completeness_;
        bool verboseForPlotting_;

        std::map<std::vector<int>, ScalarType> hashTable_;
        typename std::map<std::vector<int>, ScalarType>::iterator iter_;
};

#endif