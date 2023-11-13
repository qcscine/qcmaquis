/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MODELS_VIBRATIONAL_NMODEPAIRED_H
#define MODELS_VIBRATIONAL_NMODEPAIRED_H

#ifdef DMRG_VIBRATIONAL

#include <set>
#include <cmath>
#include <sstream>
#include <functional>
#include <numeric>
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
class NModeModelPaired : public model_impl<Matrix, TrivialGroup> {
    // Types definition
    using base = model_impl<Matrix, TrivialGroup>;
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
public:

    /**
    * @brief Class constructor
    * @param lattice object representing the DMRG lattice
    * @param parameters container with the DMRG parameters
    * @param verbose if true, prints information regarding the Hamiltonian terms
    */
    NModeModelPaired(const Lattice& lattice_, BaseParameters& parameters_, bool verbose)
            : lattice(lattice_), parameters(parameters_), tag_handler(new table_type()), phys_indexes(0) //check whether everything needed is here
    {
        // Loads in the relevant parameters
        lattice_size = parameters_["L"].as<int>();
        numModes = parameters_["nmode_num_modes"].as<int>();
        maxCouplingDegree = parameters["nmode_max_coupling"];
        //resize the vector containing the physical indices
        phys_indexes.resize(numModes);
        // Analyzes consistency of numModals parameter
        nModalsVec = parameters_["nmode_num_basis"].template as<std::vector<int> >();
        if (nModalsVec.size() == 1) {
            auto nModals = nModalsVec[0];
            nModalsVec = std::vector<int>(numModes, nModals);
        }
        else if (nModalsVec.size() != numModes) {
            throw std::runtime_error("nmode_num_basis needs to be either a single integer or a list with lenght L");
        }
        // Loads the physical indices
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        for (int iMode = 0; iMode < numModes; iMode++) {
            phys_indexes[iMode].insert(std::make_pair(C, nModalsVec[iMode]));
        }
        // Decides how many different dimensions there are
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        //comvert std::set to std::vector
        std::vector<int> nModalsUniqueVec(nModalsUnique.begin(), nModalsUnique.end());
        // Loads the vector with the site types
        siteTypes.reserve(lattice_size);
        for (int iSite = 0; iSite < lattice_size; iSite++) {
            siteTypes.push_back(lattice.get_prop<int>("type", iSite)); //this has to be compatible with the chosen lattice
        }
        // == DEFINITION OF THE ELEMENTARY OPERATORS ==
        // For each mode along with all its associated modals, 
        // we define the identity, the creation, the annihilation, and the count operator.
        // define the number of operators we need
        std::vector<op_t> ident_op, count_op, destroy_op, create_op, paired_op;
        //generate matrices
        //if all modes have the same modal basis size
        if(nModalsUnique.size() == 1){
            int overallDimension = nModalsVec[0];
            Matrix mident(overallDimension, overallDimension, 0.), mcount(overallDimension, overallDimension, 0.);
            std::vector <Matrix> mpairedVec;
            mident(0, 0) = 1.;
            for (int n = 0; n < overallDimension; n++) { //change this
                for (int m = 0; m < overallDimension; m++){
                    Matrix mpaired(overallDimension, overallDimension, 0.);
                    mpaired(n,m) = 1;
                    mpairedVec.push_back(mpaired);
                    
                }
                if(n !=0){
                    mident(n, n) = 1.;
                    mcount(n, n) = value_type(n); //same count operator as in Watson Model
                }
            }
            // Local operators
            std::vector <op_t> paired_op_locVec;
            op_t ident_op_loc, count_op_loc;
            ident_op_loc.insert_block(mident, C, C);
            count_op_loc.insert_block(mcount, C, C);
            for (int i = 0; i < mpairedVec.size(); i++){
                op_t paired_op_loc;
                paired_op_loc.insert_block(mpairedVec[i], C, C);
                paired_op_locVec.push_back(paired_op_loc);
            }
            // Updates the vectors
            ident_op.push_back(ident_op_loc);
            count_op.push_back(count_op_loc);
            for (int i = 0; i < paired_op_locVec.size(); i++){
                op_t pairedToPushBack = paired_op_locVec[i];
                paired_op.push_back(pairedToPushBack);
            }
            // Creates the final tags and update the table
            ident = modelHelper<Matrix, TrivialGroup>::register_all_types(ident_op, tag_detail::bosonic, tag_handler);
            count = modelHelper<Matrix, TrivialGroup>::register_all_types(count_op, tag_detail::bosonic, tag_handler);
            paired = modelHelper<Matrix, TrivialGroup>::register_all_types(paired_op, tag_detail::bosonic, tag_handler);
        }


        else { //if the modal bases of the various modes are of different size
            for (const auto& nModals_idx: nModalsUnique) {
                int overallDimension = nModalsVec[nModals_idx];
                Matrix mident(overallDimension, overallDimension, 0.), mcount(overallDimension, overallDimension, 0.);
                std::vector <Matrix> mpairedVec;
                mident(0, 0) = 1.;
                for (int n = 1; n < overallDimension; n++) {
                    for (int m = 1; m < overallDimension; m++){
                        Matrix mpaired(overallDimension, overallDimension, 0.);
                        mpaired(n,m) = 1; 
                        mpairedVec.push_back(mpaired);
                }
                    mident(n, n) = 1.;
                    mcount(n, n) = value_type(n); //same count operator as in Watson Model
                }
                //local operators
                std::vector <op_t> paired_op_locVec;
                op_t ident_op_loc, count_op_loc;
                ident_op_loc.insert_block(mident, C, C);
                count_op_loc.insert_block(mcount, C, C);
                for (int i = 0; i < mpairedVec.size(); i++){
                    op_t paired_op_loc;
                    paired_op_loc.insert_block(mpairedVec[i], C, C);
                    paired_op_locVec.push_back(paired_op_loc);
                }
                // Updates the vectors
                ident_op.push_back(ident_op_loc);
                count_op.push_back(count_op_loc);
                for (int i = 0; i < paired_op_locVec.size(); i++){
                    op_t pairedToPushBack = paired_op_locVec[i];
                    paired_op.push_back(pairedToPushBack);
                }
            }
            // Creates the final tags and update the table
            ident = modelHelper<Matrix, TrivialGroup>::register_all_types(ident_op, tag_detail::bosonic, tag_handler);
            count = modelHelper<Matrix, TrivialGroup>::register_all_types(count_op, tag_detail::bosonic, tag_handler);
            paired = modelHelper<Matrix, TrivialGroup>::register_all_types(paired_op, tag_detail::bosonic, tag_handler);
        }
           
    } //end of constructor

  /**
   * @brief Method to load the terms.
   * This method populates the [terms_] member with the Hamiltonian coefficients
   */

    void create_terms() override { //this is taken from nu1/model.hpp. i dont think i have to change this but check and make sure...
        std::cout << "Parsing integral file" << std::endl;
        auto HamiltonianTerms = Vibrational::detail::NModeIntegralParser<double>(parameters, lattice);
        int hamiltonianSize = HamiltonianTerms.first.size();
        std::cout << "size of vector Hamiltonian_term : " << hamiltonianSize << std::endl;
        std::cout << "Processing Second-Quantization Hamiltonian" << std::endl;
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        for (int iTerm = 0; iTerm < hamiltonianSize; iTerm++) {
            positions_type positions;
            operators_type operators;
            convertLineToOperators(HamiltonianTerms.first[iTerm], positions, operators);
            std::cout << "hello" << std::endl;
            if (positions.size()/2 <= maxCouplingDegree) {
                std::cout << "goodbye" << std::endl;
                std::cout << "positions size " << positions.size() << std::endl;
                std::cout << "operators size " << operators.size() << std::endl; 
                auto matrixElement = static_cast<value_type>(HamiltonianTerms.second[iTerm]);
                modelHelper<Matrix, TrivialGroup>::add_term(positions, operators, matrixElement, tag_handler, this->terms_);
                std::cout << "called model helper" << std::endl;
                std::cout << "position " << positions[0] << std::endl;
                std::cout << "operator " << operators[0] << std::endl;
                std::cout << "matrixElement " << matrixElement << std::endl;
            }
        }
        std::cout << "Second-Quantization Hamiltonian processed" << std::endl;
    }

    /** @brief Getter for the physical dimension of a given type */
    Index<TrivialGroup> const& phys_dim(size_t type) const override { return phys_indexes[type]; }

    /** @brief Getter for the identity operator */
    tag_type identity_matrix_tag(size_t type) const override {
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        if(nModalsUnique.size() == 0) return ident[0];
        else{
            std::set<int>::iterator it = nModalsUnique.find(nModalsVec[type]);
            if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
            int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
            return ident[indexInSet];
        }
    }

    /** @brief Getter for the count operator */
    tag_type count_matrix_tag(size_t type) const {
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        if(nModalsUnique.size() == 0) return ident[0];
        else{
            std::set<int>::iterator it = nModalsUnique.find(nModalsVec[type]);
            if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
            int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
            return count[indexInSet];
        }
    }

    /** @brief Getter for the filling operator */
    tag_type filling_matrix_tag(size_t type) const override { return identity_matrix_tag(type); }


    /** @brief Gets the quantum number associated with the wfn */
    typename TrivialGroup::charge total_quantum_numbers(BaseParameters& parms) const override {
    return typename TrivialGroup::charge();
    }

    tag_type get_operator_tag(const std::string& name, size_t type) const {
        if (name == "n")
            return count_matrix_tag(type);
        else if (name == "id")
            return identity_matrix_tag(type);
        else if (name == "fill")
            return identity_matrix_tag(type);
        //else if (name == "bdag")
            //return create[type];
        //else if (name == "b")
            //return destroy[type];
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    /** @brief Getter for the tag_handler */
    table_ptr operators_table() const override { return tag_handler; }

    /** @brief Update the model with the new parameters */
    void update(BaseParameters const &p) override {
      // TODO: update this->terms_ with the new parameters
      throw std::runtime_error("update() not yet implemented or this model.");
    }

    measurements_type measurements() const override {//TODO
        measurements_type meas;
        return meas;
    }


private:

    template<class IntegralContainer>
    void convertLineToOperators(const IntegralContainer& ham_term, positions_type& pos, operators_type& ops)
    {
        assert (ham_term.size() % 2 == 0);
        int jCont = 0;
        ops.reserve(ham_term.size());
        pos.reserve(ham_term.size());
        std::set<int> nModalsUnique(nModalsVec.begin(), nModalsVec.end());
        if (nModalsUnique.size() == 1){ //if all modes have the same size physical basis
            do {
                // Retrieves matrix element
                auto mode = ham_term[2*jCont]-1; //mode number
                auto modal = ham_term[2*jCont+1]; //modal number
                int modalToCreate;
                int modalToDestroy;
                bool pair = false;
                int basisSize = nModalsVec[0]; //without the vacuum state
                assert(mode < lattice_size);
                if (jCont % 2 == 0){
                    modalToCreate = modal; 
                }
                else{
                    modalToDestroy = modal;
                    pair = true;
                }
                if(pair){
                    int index = modalToCreate*(basisSize) + modalToDestroy;
                    ops.push_back(paired[index]);
                    pos.push_back(mode);
                    std::cout << "index" << modalToCreate*(basisSize) + modalToDestroy << std::endl;
                    std::cout << "operator : " << paired[index] << std::endl;
                    std::cout << "position in lattice : " << mode << std::endl;
                } 
                pair = false;
                jCont += 1;
            }
            while (2*jCont < ham_term.size() && ham_term[2*jCont] != -1);
        }
        else{ //if physical bases have different sizes //TODO
            do {
                auto mode = ham_term[2*jCont]-1; //mode number
                auto modal = ham_term[2*jCont+1]; //modal number
                int modalToCreate;
                int modalToDestroy;
                bool pair = false;
                assert(mode < lattice_size);
                int dimension = nModalsVec[mode];
                std::set<int>::iterator it = nModalsUnique.find(dimension); //find pointer to desired dimension in the set
                if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
                int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
                // Iterate through the set and sum up elements up to the target index
                int currentIndex = 0;
                int collector = 0;
                int element;
                for (auto it = nModalsUnique.begin(); it != nModalsUnique.end() && currentIndex < indexInSet; ++it, ++currentIndex) {
                    element = *it;
                    element *= element;
                    collector += element;
                }
                if (jCont % 2 == 0)
                    modalToCreate = modal;
                else{
                    modalToDestroy = modal; 
                    pair = true;
                }
                if(pair){
                    int index = (collector -1) modalToCreate*(basisSize) + modalToDestroy;
                    ops.push_back(paired[index]);
                    pos.push_back(mode);

                }
                pair = false;
                jCont += 1;
            }
            while (2*jCont < ham_term.size() && ham_term[2*jCont] != -1);
        }
    }





private:
    /** Ref to the lattice object */
    const Lattice& lattice;
    /** Class member indicating the highest value of the Taylor operator */
    int lattice_size, numModes, maxCouplingDegree;
    /** Parameter container */
    BaseParameters& parameters;
    /** Physical basis */
    std::vector<Index<TrivialGroup>> phys_indexes;
    /** Pointer to the tag_handler */
    std::shared_ptr<TagHandler<Matrix, TrivialGroup> >  tag_handler;
    /** Elementary Operators*/
    operators_type ident, count, create, destroy, paired;
    /**Type asscoiated with the different sites*/
    std::vector<int> siteTypes;
    /**Vector contining the number of modals per mode*/
    std::vector<int> nModalsVec;
    /**variable containing the total number of operator pairs needed per mode*/
    int numOperators;
    /**set containing the unique dimensions of the required oeprators*/
};

#endif // DMRG_VIBRATIONAL

#endif


