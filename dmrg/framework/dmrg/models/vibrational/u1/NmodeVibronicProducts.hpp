#ifdef DMRG_VIBRONIC

#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"
#include "dmrg/models/vibrational/VibronicIntegralParser.hpp"

template<class Matrix>
class ExcitonicNmode : public model_impl<Matrix, U1>
{
public:
    //Types definition
    using base = model_impl<Matrix, U1>;
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using operators_type = typename std::vector<tag_type>;
    using term_descriptor = typename base::term_descriptor;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using value_type = typename Matrix::value_type;
    using pos_t = typename Lattice::pos_t;

    /**
    * @brief Class constructor
    * @param lattice object representing the DMRG lattice
    * @param parameters container with the DMRG parameters
    */

    ExcitonicNmode(const Lattice& lat_, BaseParameters & model_)
            : lat(lat_), model(model_), L_(model["L"]), /*tag_handler(new table_type()),*/ n_ele_states_(model["vibronic_num_elestates"]), 
            n_vib_states_(model["vibronic_num_vibmodes"]), n_particles_(model["vibronic_num_molecules"]), phys_indexes(0), J_(0.),
            epsilon_(1.), only_nn_(false)
    {
        // Constructor
        tag_handler = std::make_shared<TagHandler<Matrix, U1>>();
        only_nn_ = true; // currently hardcoded
        J_ = model["vibronic_J_coupling"].as<value_type>();
        epsilon_ = model["vibronic_J_excitation"].as<value_type>();
        nMaxVec = model["Nmax"].as<std::vector<int> >();
        n_connectingmodes = model["vibronic_num_connectingmodes"].as<int>();
        num_vibtypes = n_vib_states_*n_particles_-n_connectingmodes;
        op_t  ident_ele_op, create_ele_op, destroy_ele_op, count_ele_op, count_ele_op_gs; //Electronic operators:
        
        // Analyzes consistency of nMax parameter
        if (nMaxVec.size() == 1) {
            auto nMax = nMaxVec[0];
            nMaxVec = std::vector<int>(num_vibtypes, nMax);
        }
        else if (nMaxVec.size() != num_vibtypes) {
            throw std::runtime_error("Nmax needs to be either a single integer or a list with lenght n_modes*n_particles-n_connectingmodes");
        }

        // Definition of the physical dimensions.
        // Vibrational dimensions:
        phys_indexes.resize(num_vibtypes+1); //currently only one electronic state possible
        for (int iMode = 1; iMode <= num_vibtypes; iMode++)
            phys_indexes[iMode].insert(std::make_pair(0, nMaxVec[iMode-1]));

        // Electronic dimensions:
        phys_indexes[0].insert(std::make_pair(0, 1));
        phys_indexes[0].insert(std::make_pair(1, 1));

        std::cout << "PRINTING PHYSICAL INDICES" << std::endl;
        for (const auto& iEl: phys_indexes)
            std::cout << "phys index is " << iEl << std::endl;

        // Handle electronic operators
        ident_ele_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_ele_op.insert_block(Matrix(1, 1, 1), 1, 1);
        create_ele_op.insert_block(Matrix(1, 1, 1), 0, 1); 
        destroy_ele_op.insert_block(Matrix(1, 1, 1), 1, 0);
        count_ele_op.insert_block(Matrix(1, 1, 1), 1, 1); //count for excited state
        count_ele_op_gs.insert_block(Matrix(1, 1, 1), 0, 0); //count for ground state

        std::cout << "electronic identity" << ident_ele_op << std::endl; 
        std::cout << "electronic creation" << create_ele_op << std::endl;
        std::cout << "electronic destroyer" << destroy_ele_op << std::endl;
        std::cout << "electronic count" << count_ele_op << std::endl;
        std::cout << "electronic count, ground state" << count_ele_op_gs << std::endl;

        // Register electronic operators
        ident_ele = tag_handler->register_op(ident_ele_op, tag_detail::bosonic);
        std::cout << "registered electronic identity with tag " << ident_ele << std::endl;
        create_ele = tag_handler->register_op(create_ele_op, tag_detail::bosonic);
        std::cout << "registered electronic creator with tag " << create_ele << std::endl;
        destroy_ele = tag_handler->register_op(destroy_ele_op, tag_detail::bosonic);
        std::cout << "registered electronic destroyer with tag " << destroy_ele << std::endl;
        count_ele = tag_handler->register_op(count_ele_op, tag_detail::bosonic);
        std::cout << "registered electronic count op with tag " << count_ele << std::endl;
        count_ele_gs = tag_handler-> register_op(count_ele_op_gs, tag_detail::bosonic);
        std::cout << "registered electronic ground state count op with tag " << count_ele_gs << std::endl;

        // Handle vibrational operators
        std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
        //TODO: SiteTypes?
        std::vector<op_t> ident_op, count_op, destroy_op, create_op, paired_op;
        //If modal bases have different number of modals
        for (const auto& nModals_idx: nModalsUnique) {
            int overallDimension = nModals_idx;
            std::cout << "Matrix dimension is " << overallDimension << std::endl;
            Matrix mident(overallDimension, overallDimension, 0.), mcount(overallDimension, overallDimension, 0.);
            std::vector <Matrix> mpairedVec;
            for (int n = 0; n < overallDimension; n++) {
                for (int m = 0; m < overallDimension; m++){
                    Matrix mpaired(overallDimension, overallDimension, 0.);
                    mpaired(n,m) = 1; 
                    std::cout << "created a paired operator with the entry 1 in row index " << n << " and column index " << m << std::endl;
                    mpairedVec.push_back(mpaired);
                }
                mident(n, n) = 1.;
                if (n != 0) mcount(n, n) = value_type(n); //same count operator as in Watson Model
            }
            //DEBUG
            std::cout << "printing operators" << std::endl;
            std::cout << "identity " << mident << std::endl;
            std::cout << "count " << mcount << std::endl;
            for(const auto &iEl : mpairedVec) std::cout << "paired operator " << iEl << std::endl;
            //local operators
            std::vector <op_t> paired_op_locVec;
            op_t ident_op_loc, count_op_loc;
            ident_op_loc.insert_block(mident, 0, 0);
            count_op_loc.insert_block(mcount, 0, 0);
            for (int i = 0; i < mpairedVec.size(); i++){
                op_t paired_op_loc;
                paired_op_loc.insert_block(mpairedVec[i], 0, 0);
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
        ident = modelHelper<Matrix, U1>::register_all_types(ident_op, tag_detail::bosonic, tag_handler);
        for(const auto &iEl : ident)
            std::cout << "registered a vib identity with tag " << iEl << std::endl;
        count = modelHelper<Matrix, U1>::register_all_types(count_op, tag_detail::bosonic, tag_handler);
        for(const auto &iEl : count)
            std::cout << "registered a vib count op with tag " << iEl << std::endl;
        paired = modelHelper<Matrix, U1>::register_all_types(paired_op, tag_detail::bosonic, tag_handler);
        for(const auto &iEl : paired)
            std::cout << "registered a vib paired op with tag " << iEl << std::endl;
    }

    void create_terms() override {
        std::cout << "Parsing integral file" << std::endl;
        auto hamiltonianTerms = Vibrational::detail::parseIntegralExcitonicNmode<value_type>(model, lat);
        auto hamiltonianSize = hamiltonianTerms.second.size();
        for (int i_body = 0; i_body < n_particles_; i_body++) { //loop over monomers
            std::cout << "iteration of " << i_body << "th monomer" << std::endl;
            int flag = 0;
            std::vector<int> vec_jnk(2);
            std::vector<int> vec_jnk_next(2);
            vec_jnk_next[0] = i_body+1;
            vec_jnk[0] = i_body;
            for (int idx = 0; idx < hamiltonianTerms.first.size(); idx++){ //loop over rows of integral file
                std::cout << "iteration of " << idx << "th row of integral file" << std::endl;
                int connecting = -1;
                int ele_state = -1;
                if(model["vibronic_max_coupling_nmode"].as<int>() != 3){
                    auto it = std::find(hamiltonianTerms.first[idx].begin(), hamiltonianTerms.first[idx].end(), -1);
                    if (it != hamiltonianTerms.first[idx].end()) {
                    // Calculate the index by subtracting begin() iterator from the found iterator
                        int maxIndex = std::distance(hamiltonianTerms.first[idx].begin(), it);
                        ele_state = hamiltonianTerms.first[idx][maxIndex-1];
                        std::cout << "ele_state is " << ele_state << std::endl;
                        connecting = hamiltonianTerms.first[idx][maxIndex-2];
                        std::cout << "connecting variable is " << connecting << std::endl;
                    } else {
                        throw std::runtime_error("Problem with returned operator array from integral parser");
                    }
                }
                else {
                    ele_state = hamiltonianTerms.first[idx][hamiltonianTerms.first.size()-1];
                    std::cout << "ele_state is " << ele_state << std::endl;
                    connecting = hamiltonianTerms.first[idx][hamiltonianTerms.first.size()-2];
                    std::cout << "connecting variable is " << connecting << std::endl;
                }
                if (connecting == -1 || ele_state == -1) throw std::runtime_error("Error reading variable <<connecting>> and/or <<ele_state>>");
                auto matrixElement = static_cast<value_type>(hamiltonianTerms.second[idx]);
                if(i_body == n_particles_-1 && connecting) break;
                std::vector<tag_type> operators;
                std::vector<pos_t> positions;
                std::vector<int> modes;
                std::vector<int> modals;
                bool skip = false;
                for(int i = 0; i < hamiltonianTerms.first[idx].size(); i++){ //loop over coupled modes and modals in a single integral line
                    if (hamiltonianTerms.first[idx][i] == -1 || hamiltonianTerms.first[idx][i+2] == -1) break;
                    else if (i % 2 == 0) modes.push_back(hamiltonianTerms.first[idx][i]-1); //-1 so that mode index starts at zero
                    else modals.push_back(hamiltonianTerms.first[idx][i]);
                }
                for(int i = 0; i < modes.size(); i+=2){
                    vec_jnk[1] = modes[i];
                    int localDimension = nMaxVec[i_body*n_vib_states_+modes[i]];
                    int modalToCreate = modals[i];
                    int modalToDestroy = modals[i+1];
                    // Calculating the index of the given operator in the "paired" vector
                    // Iterate through the set
                    int sumSquaredSmaller = 1;
                    std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
                    for (const auto& element : nModalsUnique) {
                        // Check if the element is smaller than the given integer
                        if (element < localDimension) {
                            int squaredElement = element * element; // Square the element
                            sumSquaredSmaller += squaredElement; // Add the squared element to the sum
                        }
                    }
                    if(modalToCreate >= localDimension || modalToDestroy >= localDimension) skip = true;
                    int indexOp = (sumSquaredSmaller-1) + modalToCreate*(localDimension) + modalToDestroy;
                    std::cout << "modalToDestroy " << modalToDestroy << std::endl;
                    std::cout << "modalToCreate " << modalToCreate << std::endl;
                    std::cout << "index " << indexOp << std::endl;
                    operators.push_back(paired[indexOp]); //define ops
                    positions.push_back(lat.get_prop<int>("vibindex", vec_jnk)); //define pos
                    std::cout << "pushed back paired operator with tag " << operators[0] << " at position " << positions[0] << std::endl; 
                }
                if(skip) continue;
                // Add electronic contribution
                // Add the count operator for the specific excited states.
                if (ele_state == 1) { //if electronic excited state potential
                    std::cout << "detected electronic excited state" << std::endl;
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(count_ele);
                    std::cout << "created |1><1| term" << std::endl;
                    //BEGIN NEW
                    if( (i_body < n_particles_-1) && connecting == 1){
                        //create |1><1||0><0| term
                        vec_jnk_next[1] = 0;
                        positions.push_back(lat.get_prop<int>("eleindex", vec_jnk_next));
                        operators.push_back(count_ele_gs);
                        modelHelper<Matrix, U1>::add_term(positions, operators, matrixElement, tag_handler, this->terms_, true);
                        for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                        for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                        std::cout << "created |1><1||0><0| term" << std::endl;
                        //create |1><1||1><1| term
                        operators.pop_back();
                        operators.push_back(count_ele);
                        modelHelper<Matrix, U1>::add_term(positions, operators, matrixElement, tag_handler, this->terms_, true);
                        for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                        for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                        std::cout << "created |1><1||1><1| term" << std::endl;
                        //create |0><0||1><1| term
                        operators.pop_back(); //remove count_ele of next site
                        operators.pop_back(); //remove count_ele of current site
                        operators.push_back(count_ele_gs);
                        operators.push_back(count_ele);
                        modelHelper<Matrix, U1>::add_term(positions, operators, matrixElement, tag_handler, this->terms_, true);
                        for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                        for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                        flag = 1;
                        std::cout << "created |0><0||1><1| term" << std::endl;
                    }
                    

                }
                else{ //if electronic ground state potential
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(count_ele_gs);
                    std::cout << "detected electronic ground state" << std::endl;
                    
                    if( (i_body < n_particles_-1) && connecting == 1){ 
                        vec_jnk_next[1] = 0;
                        positions.push_back(lat.get_prop<int>("eleindex", vec_jnk_next));
                        operators.push_back(count_ele_gs);
                        modelHelper<Matrix, U1>::add_term(positions, operators, matrixElement, tag_handler, this->terms_, true);
                        for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                        for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                        flag = 1;
                        std::cout << "mark for neighbouring excited vibration" << std::endl; 
                    }
                }
                // Builds the term of the Hamiltonian
                if( !(i_body == n_particles_-1 && connecting == 1) && flag == 0 ){ 
                    modelHelper<Matrix, U1>::add_term(positions, operators, matrixElement, tag_handler, this->terms_, true);
                    for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                    for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                }
                flag = 0;
            }
        }
        // Add J coupling
        std::vector<int> vec_jnk(2);
        for (int i1_body = 0; i1_body < n_particles_; i1_body++) {
            for (int i2_body = 0; i2_body < n_particles_; i2_body++) {
                if ((only_nn_ && (i1_body-i2_body == 1 || i2_body-i1_body == 1)) || (!only_nn_ && i1_body!=i2_body)) {
                    std::vector<tag_type> operators;
                    std::vector<pos_t> positions;
                    vec_jnk[0] = i1_body;
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    vec_jnk[0] = i2_body;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(create_ele);
                    operators.push_back(destroy_ele);
                    modelHelper<Matrix, U1>::add_term(positions, operators, J_, tag_handler, this->terms_);
                    for(const auto& iEl : positions)
                            std::cout << "added term for position " << iEl << std::endl;
                    for(const auto& iEl : operators)
                            std::cout << "with operator " << iEl << std::endl;
                }
            }
        }
        std::cout << "added J coupling" << std::endl;
    }

    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
    }

    /** @brief Getter for the physical basis */
    Index<U1> const& phys_dim(size_t type) const { return phys_indexes[type];
    }

    /** @brief Identity matrix getter */
    tag_type identity_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if(type != 0){
            std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
            std::set<int>::iterator it = nModalsUnique.find(nMaxVec[type-1]);
            if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
            int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
            return ident[indexInSet];
        }
        else
            ret = ident_ele;
        return ret;
    }

    /** @brief Count matrix getter */
    tag_type count_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if(type != 0){
            std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
            std::set<int>::iterator it = nModalsUnique.find(nMaxVec[type-1]);
            if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
            int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
            return count[indexInSet];
        }
        else
            throw std::runtime_error("No number operator for electronic sites");
        return ret;
    }

    /** @brief Filling matrix getter */
    tag_type filling_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if ((type <= num_vibtypes) && (type > 0)){
            std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
            std::set<int>::iterator it = nModalsUnique.find(nMaxVec[type-1]);
            if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
            int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
            return ident[indexInSet];
        }
        else if (type == 0)
            ret = ident_ele;
        else
            throw std::runtime_error("Site type not recognized") ;
        return ret ;
    }

    /** @brief Charge getter */
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const { return parms["vibronic_num_excitons"]; }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        //if (name == "n")
            //return count_matrix_tag(type);
        if (name == "id")
            return identity_matrix_tag(type);
        else if (name == "fill")
            return identity_matrix_tag(type);
        else
          throw std::runtime_error("Operator not valid for this model.");
        return 0;
        //creation and annihilation not yet implemented
    }

    /** @brief Getter for the operator table */
    table_ptr operators_table() const { return tag_handler; }

    measurements_type measurements() const
    {
        //Types definitions
        using op_vec = std::vector<op_t>;
        using bond_element = std::vector<std::pair<op_vec, bool> >;
        //Variable declaration
        measurements_type meas;
        std::set<int> nModalsUnique(nMaxVec.begin(), nMaxVec.end());
        //Ground State Population
        if (model.is_set("MEASURE[Population]")) {
            for (std::size_t idx = 0; idx < n_particles_; idx++){
                std::string name = "PopulationState"+std::to_string(idx);
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back((n_vib_states_+n_ele_states_)*idx);
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                for (std::size_t idx = 0; idx < num_vibtypes + n_ele_states_; idx++){
                    identities_local.push_back(this->identity_matrix(idx));
                    fillings_local.push_back(this->filling_matrix(idx));
                }
                // Bonds element (the actual operator involved in the measurement)
                bond_element ops;
                op_vec local_op_vec;
                local_op_vec.push_back(tag_handler->get_op(count_ele));
                for (std::size_t idx = 0; idx < num_vibtypes; idx++) {
                    //figure out what dimension the given site is
                    std::set<int>::iterator it = nModalsUnique.find(nMaxVec[idx]);
                    if(it == nModalsUnique.end()) std::runtime_error("Index of dimension not found in set nModalsUnique");
                    int indexInSet = std::distance(nModalsUnique.begin(), it); // extract index of entry dimension in set
                    //push back identity with correct dimensions
                    local_op_vec.push_back(tag_handler->get_op(ident[indexInSet]));
                }
                ops.push_back(std::make_pair(local_op_vec, false));
                meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local, fillings_local, ops));
            }
        }

        if (model.is_set("MEASURE[ModeExcitationDegree]")){
            int typeCount = 1;
            for(std::size_t iBody = 0; iBody < n_particles_; iBody++){ //loop over monomers
                for(std::size_t iMode = 0; iMode < n_vib_states_; iMode++){ //loop over vibrational modes
                    if( (iBody == n_particles_-1) && (iMode >= (n_vib_states_-n_connectingmodes)) ) break;
                    std::string name = "Monomer"+std::to_string(iBody)+"ExcitationMode"+std::to_string(iMode);
                    // Generates vectors for the positions
                    std::vector<pos_t> pos_internal(0);
                    std::vector<std::vector<pos_t>> pos_local(0);
                    pos_internal.push_back((iMode+1) + (iBody*(n_vib_states_+n_ele_states_)));
                    pos_local.push_back(pos_internal);
                    // Account for fillings and identities
                    op_vec identities_local, fillings_local;
                    for(std::size_t idx = 0; idx <= num_vibtypes; idx++){
                        identities_local.push_back(this->identity_matrix(idx)); 
                        fillings_local.push_back(this->filling_matrix(idx)); 
                    }
                    bond_element ops;
                    op_vec local_op_vec;
                    local_op_vec.push_back(tag_handler->get_op(ident_ele)); //electronic identity
                    for(std::size_t idx = 1; idx <= num_vibtypes; idx++){
                        if(typeCount == idx){
                            auto localOperator = tag_handler->get_op(count_matrix_tag(idx)); //TODO get correct operator
                            local_op_vec.push_back(localOperator);
                        }
                        else local_op_vec.push_back(tag_handler->get_op(identity_matrix_tag(idx))); 
                    }
                    typeCount++;
                    ops.push_back(std::make_pair(local_op_vec, false));
                    meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local, fillings_local, ops));
                }
            }
        }
        return meas;
    }


private:
    const Lattice& lat;
    BaseParameters& model;
    value_type J_, epsilon_;
    bool only_nn_;
    std::size_t L_;
    std::size_t n_ele_states_, n_vib_states_, n_particles_;
    std::vector<Index<U1>> phys_indexes;
    std::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    std::vector<int> nMaxVec;
    int n_connectingmodes;
    int num_vibtypes;
    //operators
    operators_type ident, count, create, destroy, paired; //vibrational operators
    tag_type ident_ele, count_ele, count_ele_gs, create_ele, destroy_ele; //electronic operators


    


};

#endif // DMRG_VIBRONIC
