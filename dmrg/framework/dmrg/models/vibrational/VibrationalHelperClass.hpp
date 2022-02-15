/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2021 Laboratory for Physical Chemistry, ETH Zurich
 *               2021 by Alberto Baiardi <abaiardi@phys.ethz.ch>
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

#ifndef VIBRATIONAL_MODEL_HELPER_CLASS_H
#define VIBRATIONAL_MODEL_HELPER_CLASS_H

template<class Matrix, class SymmGroup>
class VibrationalHelpers {
public:
    // TYpes definition
    using OperatorType = typename operator_selector<Matrix, SymmGroup>::type;
    
    /** @brief Generates q^n terms */
    static std::vector<OperatorType> generatePowersOfPositionOperator(int maxCoupling, int nMax, const OperatorType& identityOp,
                                                                      const OperatorType& positionOp)
    {
        // -- Creates the powers of the position/momentum operator --
        std::vector<OperatorType> powersOfPositions(maxCoupling+1);
        powersOfPositions[0] = identityOp;
        powersOfPositions[0].resize_block(0, nMax, nMax);
        OperatorType q = identityOp;
        for (int iOrder = 0; iOrder < maxCoupling; iOrder++) 
        {
            OperatorType tmpQ;
            gemm(q, positionOp, tmpQ);
            powersOfPositions[iOrder+1] = tmpQ;
            q = tmpQ;
            powersOfPositions[iOrder+1].resize_block(0, nMax, nMax);
        }
        return powersOfPositions;
    }

    /** @brief Generates p^n terms */
    static std::vector<OperatorType> generatePowersOfMomentumOperator(int maxCoupling, int nMax, const OperatorType& identityOp,
                                                                      const OperatorType& momentumOp)
    {
        // -- Creates the powers of the position/momentum operator --
        std::vector<OperatorType> powersOfMomentum(maxCoupling+1);
        powersOfMomentum.resize(maxCoupling+1);
        powersOfMomentum[0] = identityOp;
        powersOfMomentum[0].resize_block(0, nMax, nMax);
        OperatorType p = identityOp;
        for (int iOrder = 0; iOrder < maxCoupling; iOrder++) 
        {
            OperatorType tmpP;
            gemm(p, momentumOp, tmpP);
            powersOfMomentum[iOrder+1] = tmpP;
            p = tmpP;
            powersOfMomentum[iOrder+1].resize_block(0, nMax, nMax);
        }
        return powersOfMomentum;
    }

}; // VibrationalHelpers

#endif // VIBRATIONAL_MODEL_HELPER_CLASS_H