/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef VIBRATIONAL_MODEL_HELPER_CLASS_H
#define VIBRATIONAL_MODEL_HELPER_CLASS_H

/** @brief Enum class representing the coordinate type */
enum class WatsonCoordinateType { CartesianNormalModes, InternalNormalModes };

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
