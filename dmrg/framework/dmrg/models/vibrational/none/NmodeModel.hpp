/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MODELS_VIBRATIONAL_NONE_H
#define MODELS_VIBRATIONAL_NONE_H

#ifdef DMRG_VIBRATIONAL

#include <set>
#include <sstream>
#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalIntegralParser.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"

/**
 * @brief Class implementing the nmode vibrational Hamiltonian.
 *
 * In this model, we use the nmode quantization to map the Born-Oppenheimer
 * vibrational Hamiltonian onto the DMRG lattice.
 * In this case, all modals of a unique mode are mapped to a single site of the dmrg lattice, 
 * involving multiple creation/annihilation operator pairs per site
 */

template<class Matrix>
class NModeModel : public model_impl<Matrix, TrivialGroup> {
    // Types definition
    using base = model_impl<Matrix, NU1>;
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using terms_type = typename std::vector<term_descriptor>;
    using op_t = typename base::op_t;
    using operators_type = typename std::vector<tag_type>;
    using measurements_type = typename base::measurements_type;
    using pos_t = typename Lattice::pos_t;
    using positions_type = typename std::vector<pos_t>;
    using value_type = typename Matrix::value_type;
    using charge_type = typename NU1::charge;
public:

    /**
    * @brief Class constructor
    * @param lattice object representing the DMRG lattice
    * @param parameters container with the DMRG parameters
    * @param verbose if true, prints information regarding the Hamiltonian terms
    */
    NModeModel(const Lattice& lattice_, BaseParameters& parameters_, bool verbose)
            : lattice(lattice_), parameters(parameters_), tag_handler(new table_type()), physIndices_(0) //check whether everything needed is here
  {
        // Loads in the relevant parameters
        numModes, lattice_size = parameters["L"];
        maxCouplingDegree = parameters["nmode_max_coupling"];
        physIndices_.resize(numModes);
        //Add here the equivalent of charge types in model.hpp -> not even sure we need this
        // Analyzes consistency of numModals parameter
        numModalsVec = parameters_["nmode_num_basis"].template as<std::vector<int> >();
        if (numModalsVec.size() == 1) {
            auto nModals = numModalsVec[0];
            numModalsVec = std::vector<int>(numModes, nModals);
        }
        else if (numModalsVec.size() != numModes) {
            throw std::runtime_error("Nmax needs to be either a single integer or a list with lenght L");
        }
        // Loads the physical indices
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        for (int iMode = 0; iMode < numModes; iMode++)
        physIndices_[iMode].insert(std::make_pair(C, numModalsVec[iMode]));
        // Decides how many different dimensions there are
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        // Loads the vector with the site types
        siteTypes.reserve(lattice_size);
        for (int iSite = 0; iSite < lattice_size; iSite++)
            siteTypes.push_back(lattice.get_prop<int>("type", iSite));
        // == DEFINITION OF THE ELEMENTARY OPERATORS ==
        // For each mode along with all its associated modals, 
        // we define the identity, the creation, the annihilation, and the count operator.
        // define the number of operators we need
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        if (numModalsVec.size() == 1) auto numOperators = numModes*numModalsVec[0];
        else auto numOperators = std::reduce(numModalsVec.begin(), numModalsVec.end());
        std::vector<op_t> ident_op, create_op, destroy_op, count_op;
        ident.reserve(numOperators);
        create.reserve(numOperators);
        destroy.reserve(numOperators);
        count.reserve(numOperators);
        //gerenate matrices
        //if all modes have the same modal basis size
        if(numModalsVec.size() == 1){
            int overallDimension = numModalsVec[0];
            Matrix mcreate(overallDimension, overallDimension, 0.), mdestroy(overallDimension, overallDimension, 0.);
            Matrix mident(overallDimension, overallDimension, 0.), mcount(overallDimension, overallDimension, 0.);
            mident(0, 0) = 1.;
            for (int n = 1; n < overallDimension; n++) {
                mcreate(n-1, n) = std::sqrt(value_type(n));
                mdestroy(n-1, n) = std::sqrt(value_type(n));
                mident(n, n) = 1.;
                mcount(n, n) = value_type(n);
            }
            for (int idx = 0; idx < numOperators; idx++) {
              // Local operators
              op_t ident_op_loc, create_op_loc, destroy_op_loc, count_op_loc;
              ident_op_loc.insert_block(mident, C, C);
              create_op_loc.insert_block(mcreate, C, C);
              destroy_op_loc.insert_block(mdestroy, C, C);
              count_op_loc.insert_block(mcount, C, C);
              // Updates the vectors
              ident_op.push_back(ident_op_loc);
              create_op.push_back(create_op_loc);
              count_op.push_back(count_op_loc);
              destroy_op.push_back(destroy_op_loc);
            }
            // Creates the final tags and update the table
            ident = modelHelper<Matrix, TrivialGroup>::register_all_types(ident_op, tag_detail::bosonic, tag_handler);
            create = modelHelper<Matrix, TrivialGroup>::register_all_types(create_op, tag_detail::bosonic, tag_handler);
            destroy = modelHelper<Matrix, TrivialGroup>::register_all_types(destroy_op, tag_detail::bosonic, tag_handler);
            count = modelHelper<Matrix, TrivialGroup>::register_all_types(count_op, tag_detail::bosonic, tag_handler);
            // Registers the hermitian pairs
            modelHelper<Matrix, NU1>::registerHermitianConjugates(create, destroy, tag_handler);
        }
        else { //if the modal bases of the various modes are of different size
            for (const auto& nModals_idx: nModalsUnique) {
                int overallDimension = numModalsVec[nModals_idx];
                Matrix mcreate(overallDimension, overallDimension, 0.), mdestroy(overallDimension, overallDimension, 0.);
                Matrix mident(overallDimension, overallDimension, 0.), mcount(overallDimension, overallDimension, 0.);
                mident(0, 0) = 1.;
                for (int n = 1; n < overallDimension; n++) {
                    mcreate(n-1, n) = std::sqrt(value_type(n));
                    mdestroy(n-1, n) = std::sqrt(value_type(n));
                    mident(n, n) = 1.;
                    mcount(n, n) = value_type(n);
                }
                for (int idx = 0; idx < numOperators; idx++) {
                    // Local operators
                    op_t ident_op_loc, create_op_loc, destroy_op_loc, count_op_loc;
                    ident_op_loc.insert_block(mident, C, C);
                    create_op_loc.insert_block(mcreate, C, C);
                    destroy_op_loc.insert_block(mdestroy, C, C);
                    count_op_loc.insert_block(mcount, C, C);
                    // Updates the vectors
                    ident_op.push_back(ident_op_loc);
                    create_op.push_back(create_op_loc);
                    count_op.push_back(count_op_loc);
                    destroy_op.push_back(destroy_op_loc);
                }
            }
            // Creates the final tags and update the table
            ident = modelHelper<Matrix, TrivialGroup>::register_all_types(ident_op, tag_detail::bosonic, tag_handler);
            create = modelHelper<Matrix, TrivialGroup>::register_all_types(create_op, tag_detail::bosonic, tag_handler);
            destroy = modelHelper<Matrix, TrivialGroup>::register_all_types(destroy_op, tag_detail::bosonic, tag_handler);
            count = modelHelper<Matrix, TrivialGroup>::register_all_types(count_op, tag_detail::bosonic, tag_handler);
            // Registers the hermitian pairs
            modelHelper<Matrix, NU1>::registerHermitianConjugates(create, destroy, tag_handler);
        }
  }



private:
    const Lattice& lattice;
    int lattice_size, num_modes, maxCouplingDegree;
    BaseParameters& parameters;
    std::vector<Index<NU1> > phys_indexes;
    std::shared_ptr<TagHandler<Matrix, NU1> >  tag_handler;
    operators_type ident, create, destroy, count;
    std::vector<int> siteTypes;
    std::vector<int> numModalsVec;
};

#endif // DMRG_VIBRATIONAL

#endif


