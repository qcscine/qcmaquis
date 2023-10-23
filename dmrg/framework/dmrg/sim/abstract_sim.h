/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */
#ifndef ABSTRACT_SIM_H
#define ABSTRACT_SIM_H

#include <map>
#include <vector>
#include "integral_interface.h"
#include "dmrg/utils/results_collector.h"

class abstract_sim {
public:
    virtual ~abstract_sim() {}
    virtual void run(const std::string& runType) = 0;
};

/**
 * @brief Virtual class defining the scheleton of a simulation object.
 * 
 * This virtual class defines the skeleton of a simulation object 
 * that runs a DMRG calculation.
 * 
 * @tparam Matrix class defining the types of the matrix entering the definition
 * of the block_matrix.
 */

template <class Matrix>
class abstract_interface_sim {
public:
    // warning, these types are defiled in model_impl already
    using meas_with_results_type = std::pair<std::vector<std::vector<int> >, std::vector<typename Matrix::value_type> >;
    using results_map_type = std::map<std::string, meas_with_results_type>;
    using RealType = typename maquis::traits::real_type<Matrix>::type;

    virtual ~abstract_interface_sim() {}
    virtual void run(const std::string& runType) = 0;
    virtual void run_measure() = 0;
    virtual RealType get_energy() = 0;
    virtual RealType getFEASTEnergy(int iState) const = 0;
    virtual results_collector& get_iteration_results() = 0;
    virtual int get_last_sweep() = 0;
    virtual results_map_type measure_out() =0;
    virtual void update_integrals(const chem::integral_map<typename Matrix::value_type> &)=0;
    virtual typename Matrix::value_type get_overlap(const std::string &) = 0;
    virtual typename Matrix::value_type getCICoefficient(std::string ciVector) = 0;
//  virtual std::string ... get_fiedler_order
};

#endif
