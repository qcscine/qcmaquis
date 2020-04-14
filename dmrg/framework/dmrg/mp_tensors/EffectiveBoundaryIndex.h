/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2020 Institute for Theoretical Physics, ETH Zurich
 *               2020 by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef EFFECTIVE_BOUNDARY_INDEX_H
#define EFFECTIVE_BOUNDARY_INDEX_H

#include <mpi_interface.h>
#include <memory>

/**
 * @brief Class containing the index of the boundaries "owned" by a given thread.
 */

class EffectiveBoundaryIndex {
public:
    // Types declaration
    //using InputType = std::unique_ptr<Scine::Mpi::MpiInterface>;
    using InputType = Scine::Mpi::MpiInterface*;

    EffectiveBoundaryIndex() : hasMPI_(false), numberOfThreads_(0) {};

    /**
     * @brief Constructor that takes as input just the 
     * @param idxMax 
     */
    EffectiveBoundaryIndex(InputType inputInterface) : hasMPI_(false), numberOfThreads_(0) {
        this->initialize(inputInterface);
    }

    /**
     * @brief Gets the indexes between 0 and idMax associated to a given thread.
     * @param idxMax upper bound for the index.
     * @return std::vector<int> vector with the indexes that "belong" to a given thread.
     */
    std::vector<int> getEffectiveIndex(int idxMax) const {
        std::vector<int> effectiveIndexes;
        int iMin=0;
        int iMax=idxMax-1;
        int remainder = idxMax&numberOfThreads_;
        bool skip = false;
        if (hasMPI_ && numberOfThreads_ > 1) {
            if (idxMax > numberOfThreads_) {
                if (remainder != 0) {
                    int integerDivision = idxMax/(numberOfThreads_-1);
                    iMin = integerDivision*idxOfThread_;
                    if (idxOfThread_ != numberOfThreads_-1) {
                        iMax = integerDivision*(idxOfThread_+1);
                    }
                    else {
                        iMax = idxMax-1;
                    }
                }
                else {
                    int integerDivision = idxMax/numberOfThreads_;
                    iMin = idxOfThread_*integerDivision;
                    iMax = (idxOfThread_+1)*integerDivision;
                }
            }
            else {
                if (idxOfThread_ < idxMax) {
                    iMin = idxOfThread_;
                    iMax = idxOfThread_+1;
                }
                else {
                    skip = true;
                }
            }
        }
        // Fills the vector
        if (!skip) {
            effectiveIndexes.reserve(iMax-iMin);
            for (int iIndex = iMin; iIndex < iMax; iIndex++)
                effectiveIndexes.push_back(iIndex);
        }
        return effectiveIndexes;
    }

    /**
     * @brief Initialize the parameter for a given MPI interface
     * @param inputInterface raw pointer to an MPI interface object.
     */
    void initialize(InputType inputInterface) {
        if ((*inputInterface).isMPIAvailable()) {
            hasMPI_ = true;
            numberOfThreads_ = (*inputInterface).getGlobalCommunicatorSize();
            idxOfThread_ = (*inputInterface).getGlobalRank();
        }
    }

private:
    /* Flag keeping track if MPI management is active or not */
    bool hasMPI_;
    /* Overall number of threads in MPI (0 if MPI is not activated) */
    int numberOfThreads_;
    /* Index of the thread that created the object */
    int idxOfThread_;
};

#endif
