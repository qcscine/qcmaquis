#ifdef DMRG_VIBRONIC

#include "dmrg/models/model_helper.hpp"
#include "dmrg/models/vibrational/VibrationalHelperClass.hpp"
#include "dmrg/models/vibrational/VibronicIntegralParser.hpp"

template<class Matrix>
class HolsteinbHubbardExcitonicExtendedHamiltonian : public model_impl<Matrix, U1>
{
public:
    //Types definition
    using base = model_impl<Matrix, U1>;
    using table_type = typename base::table_type;
    using table_ptr = typename base::table_ptr;
    using tag_type = typename base::tag_type;
    using term_descriptor = typename base::term_descriptor;
    using op_t = typename base::op_t;
    using measurements_type = typename base::measurements_type;
    using value_type = typename Matrix::value_type;
    using pos_t = typename Lattice::pos_t;

    /** 
     * @brief Model representing an extended Holstein-Hubbard Hamiltonian
     * In contrast to the Regular Holstein Hubbard Hamiltonian, excited states can be equipped
     *  with an energy offsest and non-local potentials can be present
     * Moreover, off-diagonal coordinate-independent electronic coupling terms are present.
     */
    HolsteinbHubbardExcitonicExtendedHamiltonian (const Lattice& lat_, BaseParameters & model_) 
        : lat(lat_), model(model_), tag_handler(new table_type()), L_(model["L"]), n_ele_states_(model["vibronic_num_elestates"]),
          n_vib_states_(model["vibronic_num_vibmodes"]), n_particles_(model["vibronic_num_molecules"]), phys_indexes(0), J_(0.),
          epsilon_(1.), only_nn_(false)
    {
        //constructor -- to be written
        maxCoupling = 2; // hardcoded for the moment -> should be made dynamic (TODO)
        only_nn_ = true; // hardcoded for the moment. not sure what to do with that (TODO)
        J_ = model["vibronic_J_coupling"].as<value_type>();
        epsilon_ = model["vibronic_J_excitation"].as<value_type>();
        nMaxVec = model["Nmax"].as<std::vector<int> >();
        n_connectingmodes = model["vibronic_num_connectingmodes"].as<int>();
        num_vibtypes = n_vib_states_*n_particles_-n_connectingmodes;
        check_maxVibMode = 0;
        //  Note that, since we have two different types of sites, we also have different sets of
        //  operators for electrons and for the vibrations
        op_t ident_ele_op;
        op_t create_ele_op, destroy_ele_op, count_ele_op, count_ele_op_gs;
        // Analyzes consistency of nMax parameter
        if (nMaxVec.size() == 1) {
            auto nMax = nMaxVec[0];
            nMaxVec = std::vector<int>(num_vibtypes, nMax);
        }
        else if (nMaxVec.size() != num_vibtypes) {
            throw std::runtime_error("Nmax needs to be either a single integer or a list with lenght n_modes*n_particles-n_connectingmodes");
        }
        for (int i = 0; i < nMaxVec.size(); i++){
            std::cout << "DEBUG_NMAX :" << nMaxVec[i] << std::endl;
        }
        // Definition of the physical dimensions.
        // First we manage the dimensions for the vibrations
        phys_indexes.resize(num_vibtypes+1); //TODO: Number of excitons hardcoded to 1 atm
        for (int iMode = 1; iMode <= num_vibtypes; iMode++)
            phys_indexes[iMode].insert(std::make_pair(0, nMaxVec[iMode-1]));
        //manage electronic dimensions
        phys_indexes[0].insert(std::make_pair(0, 1));
        phys_indexes[0].insert(std::make_pair(1, 1));
        // Registering the electronic operators
        ident_ele_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_ele_op.insert_block(Matrix(1, 1, 1), 1, 1);
        create_ele_op.insert_block(Matrix(1, 1, 1), 0, 1); //contrary to usual matrix index notation, the first index here corresponds to the column and the second one to the row
        destroy_ele_op.insert_block(Matrix(1, 1, 1), 1, 0); //contrary to usual matrix index notation, the first index here corresponds to the column and the second one to the row
        count_ele_op.insert_block(Matrix(1, 1, 1), 1, 1);
        count_ele_op_gs.insert_block(Matrix(1, 1, 1), 0, 0); //electronic ground state count operator
        // Creation of operator tag table for the electronic operators
        ident_ele = tag_handler->register_op(ident_ele_op, tag_detail::bosonic);
        create_ele = tag_handler->register_op(create_ele_op, tag_detail::bosonic);
        destroy_ele = tag_handler->register_op(destroy_ele_op, tag_detail::bosonic);
        count_ele = tag_handler->register_op(count_ele_op, tag_detail::bosonic);
        count_ele_gs = tag_handler -> register_op(count_ele_op_gs, tag_detail::bosonic);
        // Registering the vibrational operators
        // Decides how many different dimensions there are
        std::set<int> nMaxUnique(nMaxVec.begin(), nMaxVec.end());
        // Loop over all modes
        //NEW
        for (const auto& nMax: nMaxUnique) {
            op_t ident_vib_op, position_vib_op, momentum_vib_op;
            int overallDimension = nMax + maxCoupling;
            Matrix mpos(overallDimension, overallDimension, 0.), mmom(overallDimension, overallDimension, 0.), mident(overallDimension, overallDimension, 0.); //dimension right?
            //loads matrices
            mident(0,0) = 1.;
            for (int n=1; n < overallDimension; ++n) {
                mident(n,n) = 1.;
                mpos(n-1,n) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
                mpos(n,n-1) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
                mmom(n-1,n) = std::sqrt(value_type(n))/std::sqrt(value_type(2.));
                mmom(n,n-1) = -std::sqrt(value_type(n))/std::sqrt(value_type(2.));
            }
            position_vib_op.insert_block(mpos, 0, 0);
            momentum_vib_op.insert_block(mmom, 0, 0);
            ident_vib_op.insert_block(mident, 0, 0);
            auto powersOfPositions_op = VibrationalHelpers<Matrix, U1>::generatePowersOfPositionOperator(maxCoupling, nMax, ident_vib_op, position_vib_op);
            auto powersOfMomentum_op = VibrationalHelpers<Matrix, U1>::generatePowersOfMomentumOperator(maxCoupling, nMax, ident_vib_op, momentum_vib_op);
            ident_vib_op.resize_block(0, nMax, nMax);
            ident_vib[nMax] = tag_handler->checked_register(ident_vib_op, tag_detail::bosonic);
            std::cout << ident_vib_op << std::endl;
            std::cout << "Identity registered with tag " << ident_vib[nMax].first << " and coefficient " << ident_vib[nMax].second << std::endl;
            positionPowers[nMax].resize(maxCoupling+1);
            momentumPowers[nMax].resize(maxCoupling+1);
            positionPowers[nMax][0] = ident_vib[nMax];
            momentumPowers[nMax][0] = ident_vib[nMax];
            for (int iOrder = 1; iOrder <= maxCoupling; iOrder++) {
                std::cout << powersOfPositions_op[iOrder] << std::endl;
                std::cout << powersOfMomentum_op[iOrder] << std::endl;
                //if (powersOfPositions_op[iOrder].norm() > 1.0E-16)
                    positionPowers[nMax][iOrder] = tag_handler->checked_register(powersOfPositions_op[iOrder], tag_detail::bosonic);
                //else 
                //    positionPowers[nMax][iOrder] = ident_vib[nMax], value_type(1.);
                std::cout << "Position registered with tag " << positionPowers[nMax][iOrder].first << " and coeff " << positionPowers[nMax][iOrder].second << std::endl;
                //if (powersOfMomentum_op[iOrder].norm() > 1.0E-16)
                    momentumPowers[nMax][iOrder] = tag_handler->checked_register(powersOfMomentum_op[iOrder], tag_detail::bosonic);
                //else 
                //    momentumPowers[nMax][iOrder] = ident_vib[nMax], value_type(1.);
                std::cout << "Momentum registered with tag " << momentumPowers[nMax][iOrder].first << " and coeff " << momentumPowers[nMax][iOrder].second << std::endl;
            }   

        }
        
    }
    
    void create_terms() override { 
        // == Definition of the Hamiltonian ==
        auto hamiltonianTerms = Vibrational::detail::parseIntegralExcitonicExtended<value_type>(model, lat);
        // == Extracts whether we want to freeze any degrees of freedom
        // == Main loop over the monomers ==
        // We first loop over the number of molecules of the aggregate, and then over the
        // terms entering the vibronic Hamiltonian.
        for (int i_body = 0; i_body < n_particles_; i_body++) { //loop over all monomers
            std::vector<int> vec_jnk(maxCoupling);
            std::vector<int> vec_jnk_next(maxCoupling); //NEW
            vec_jnk[0] = i_body;
            vec_jnk_next[0] = i_body+1; //NEW
            int flag = 0; //NEW
            for (int idx = 0; idx < hamiltonianTerms.first.size(); idx++){ //loop over all rows of the integral file
                if(i_body == n_particles_-1 && abs(hamiltonianTerms.first[idx][2]) > n_vib_states_-n_connectingmodes) break;
                // Prepares the vectors to be employed when building the Hamiltonian
                std::vector<tag_type> operators;
                std::vector<pos_t> positions;
                int ele_state = hamiltonianTerms.first[idx][0]; //store wheter the parameters read in correspond to an excited or ground electronic state
                int connecting = hamiltonianTerms.first[idx][1]; //stores information wether the vibrational mode is monomer-internal or connecting two monomers
                int mode = abs(hamiltonianTerms.first[idx][2])-1;
                if (mode > check_maxVibMode) check_maxVibMode = mode+1; 
                std::cout << "MaxVibMode " << check_maxVibMode << std::endl;
                auto scalingFactor = hamiltonianTerms.second[idx];
                //Add vibrational contribution
                std::vector<int> tmpVec;
                for (int op_vib = 2; op_vib < maxCoupling+2 ; op_vib++)
                    tmpVec.push_back(hamiltonianTerms.first[idx][op_vib]);
                std::set<int> uniqueIndices = std::set<int>(tmpVec.begin(), tmpVec.end());
                //for (int op_vib = 2; op_vib < maxCoupling+2 ; op_vib++){ //loop over all operator entries in a single row of the integral file. starts at op_vib=2, bc first two entries do not encode operators.
                for (const int& index: uniqueIndices) {
                    //if (hamiltonianTerms.first[idx][op_vib] < 0){ //if momentum operator
                    int countOccurrences = std::count(tmpVec.begin(), tmpVec.end(), index);
                    if (index < 0){ //if momentum operator
                        operators.push_back(momentumPowers[nMaxVec[i_body*n_vib_states_+mode]][countOccurrences].first);
                        std::cout << "Scaling factor before " << scalingFactor << std::endl;
                        scalingFactor *= momentumPowers[nMaxVec[i_body*n_vib_states_+mode]][countOccurrences].second;
                        std::cout << "Scaling factor after " << scalingFactor << std::endl;
                        //vec_jnk[1] = -hamiltonianTerms.first[idx][op_vib]-1;
                        vec_jnk[1] = -index-1;
                        positions.push_back(lat.get_prop<int>("vibindex", vec_jnk));
                    }
                    else if (index > 0){ //if position operator
                        operators.push_back(positionPowers[nMaxVec[i_body*n_vib_states_+mode]][countOccurrences].first);
                        std::cout << "Scaling factor before " << scalingFactor << std::endl;
                        scalingFactor *= positionPowers[nMaxVec[i_body*n_vib_states_+mode]][countOccurrences].second;
                        std::cout << "Scaling factor after " << scalingFactor << std::endl;
                        //vec_jnk[1] = hamiltonianTerms.first[idx][op_vib]-1;
                        vec_jnk[1] = index-1;
                        positions.push_back(lat.get_prop<int>("vibindex", vec_jnk));  
                    }
                }
                // Add electronic contribution
                // Add the count operator for the specific excited states.
                if (ele_state == 1) { //if electronic excited state potential
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(count_ele);
                    //BEGIN NEW
                    if( (i_body < n_particles_-1) && connecting == 1){
                        //create |1><1||0><0| term
                        vec_jnk_next[1] = 0;
                        positions.push_back(lat.get_prop<int>("eleindex", vec_jnk_next));
                        operators.push_back(count_ele_gs);
                        modelHelper<Matrix, U1>::add_term(positions, operators, scalingFactor, tag_handler, this->terms_, true);
                        //create |1><1||1><1| term
                        operators.pop_back();
                        operators.push_back(count_ele);
                        modelHelper<Matrix, U1>::add_term(positions, operators, scalingFactor, tag_handler, this->terms_, true);
                        //create |0><0||1><1| term
                        operators.pop_back(); //remove count_ele of next site
                        operators.pop_back(); //remove count_ele of current site
                        operators.push_back(count_ele_gs);
                        operators.push_back(count_ele);
                        modelHelper<Matrix, U1>::add_term(positions, operators, scalingFactor, tag_handler, this->terms_, true);
                        maquis::cout << "DEBUG: created term for monomer " << i_body << " and integral file line " << idx << std::endl;
                        maquis::cout << "DEBUG: entered if statement" << std::endl;
                        maquis::cout << "DEBUG: positions: " << vec_jnk[0] << " " << vec_jnk_next[0] << " " << vec_jnk[1] << " " << vec_jnk_next[1] << std::endl;
                        flag = 1;
                    }
                    //END NEW

                }
                else{ //if electronic ground state potential
                    vec_jnk[1] = 0;
                    positions.push_back(lat.get_prop<int>("eleindex", vec_jnk));
                    operators.push_back(count_ele_gs);
                    //BEGIN NEW
                    
                    if( (i_body < n_particles_-1) && connecting == 1){ 
                        vec_jnk_next[1] = 0;
                        positions.push_back(lat.get_prop<int>("eleindex", vec_jnk_next));
                        operators.push_back(count_ele_gs);
                        modelHelper<Matrix, U1>::add_term(positions, operators, scalingFactor, tag_handler, this->terms_, true);
                        flag = 1;
                        maquis::cout << "DEBUG: created term for monomer " << i_body << " and integral file line " << idx << std::endl;
                        maquis::cout << "DEBUG: entered if statement" << std::endl; 
                    }
                    
                    //END NEW
                }
                // Builds the term of the Hamiltonian
                if( !(i_body == n_particles_-1 && connecting == 1) && flag == 0 ){ //is this check correct?
                    modelHelper<Matrix, U1>::add_term(positions, operators, scalingFactor, tag_handler, this->terms_, true);
                    maquis::cout << "DEBUG: created term for monomer " << i_body << " and integral file line " << idx << std::endl;
                }
                flag = 0;
            }                
            if(check_maxVibMode < n_vib_states_) throw std::runtime_error("more vibronic_num_vibmodes than modes in FCIDUMP file");
        }

           // Add the J term to the Hamiltonian
        std::vector<int> vec_jnk(2);
        for (int i1_body = 0; i1_body < n_particles_; i1_body++) {
            for (int i2_body = 0; i2_body < n_particles_; i2_body++) {
                if (only_nn_ && (i1_body-i2_body == 1 || i2_body-i1_body == 1) || !only_nn_ && i1_body!=i2_body) {
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
                }
            }
        }
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
        if(type == 0)
            ret = ident_ele;
        else
            ret = ident_vib.at(nMaxVec[type-1]).first;
        return ret;
    }

    /** @brief Filling matrix getter */
    tag_type filling_matrix_tag(size_t type) const
    {
        tag_type ret ;
        if ((type <= num_vibtypes) && (type > 0))
          ret = ident_vib.at(nMaxVec[type-1]).first;
        else if (type == 0)
          ret = ident_ele;
        else
          throw std::runtime_error("Site type not recognized") ;
        return ret ;
    }

    /** @brief Identity matrix getter */
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const { return parms["vibronic_num_excitons"]; }

    /** @brief Getter for the operator associated with a given string */
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "q")
            return positionPowers.at(nMaxVec[type-1])[1].first;
        else if (name == "p")
            return momentumPowers.at(nMaxVec[type-1])[1].first;
        else if (name == "aplus")
            return create_ele;
        else if (name == "a")
            return destroy_ele;
        else if (name == "id")
            return identity_matrix_tag(type);
        else if (name == "fill")
            return identity_matrix_tag(type);
        else
          throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    /** @brief Getter for the operator table */
    table_ptr operators_table() const { return tag_handler; }

    // Possible measurements associated to the model
    measurements_type measurements() const
    {
        // Types definition
        using op_vec = std::vector<op_t>;
        using bond_element = std::vector<std::pair<op_vec, bool> >;
        // Variable declaration
        measurements_type meas;
        std::set<int> nMaxUnique(nMaxVec.begin(), nMaxVec.end());
        // Ground state population
        if (model.is_set("MEASURE[Population]")) {
            for (std::size_t idx = 0; idx < n_particles_; idx++) {
                std::string name = "PopulationState"+std::to_string(idx);
                // Generates vectors for the position operators
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back((n_vib_states_+n_ele_states_)*idx); 
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                for (std::size_t idx1 = 0; idx1 <= num_vibtypes; idx1++) {
                    identities_local.push_back(this->identity_matrix(idx1));
                    fillings_local.push_back(this->filling_matrix(idx1));
                }
                // Bonds element (the actual operator involved in the measurement)
                bond_element ops;
                op_vec local_op_vec;
                local_op_vec.push_back(tag_handler->get_op(count_ele));
                for (std::size_t idx1 = 0; idx1 < num_vibtypes; idx1++) {
                    local_op_vec.push_back(tag_handler->get_op(ident_vib.at(nMaxVec[idx1]).first));
                }
                ops.push_back(std::make_pair(local_op_vec, false));
                meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local,
                                                                      fillings_local, ops));
            }
        } 
        
    if(model.is_set("MEASURE[Displacement]")){
        int n_connectingmodes = model["vibronic_num_connectingmodes"].as<int>();
        int typeCount = 0; //typecount = type - 1 (for vibrational sites) 
        for (std::size_t idx = 0; idx < n_particles_; idx++){ //loop over monomers
            for(std::size_t idx1 = 0; idx1 < n_vib_states_; idx1++){ //loop over vibrational modes
            //for(std::size_t idx1 = 0; idx_1 < (n_ele_states_+n_vib_states)) //DEBUG
                //if non existing lattice site: break. 
                //IMPORTANT: assumes underlying lattice sorting -> may be problematic
                //Assumed lattice sorting: intertwined with all connecting modes active. 
                //After an electronic site first come the local modes followed by the connecting modes
                if( (idx == n_particles_-1) && (idx1 >= (n_vib_states_-n_connectingmodes)) ) break; 
                std::string name = "Displacement"+std::to_string(idx)+"Mode"+std::to_string(idx1);
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back((idx1+1) + (idx*(n_vib_states_+n_ele_states_))); //pushes back positions of vibrational states
                maquis::cout << "pos of pushed back positions: " << (idx1+1) + (idx*(n_vib_states_+n_ele_states_)) << std::endl;
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                for (std::size_t idx2 = 0; idx2 <= num_vibtypes; idx2++) {
                    identities_local.push_back(this->identity_matrix(idx2));
                    fillings_local.push_back(this->filling_matrix(idx2));
                }
                bond_element ops;
                op_vec local_op_vec;
                local_op_vec.push_back(tag_handler->get_op(ident_ele));
                for (std::size_t idx2 = 0; idx2 < num_vibtypes; idx2++) {
                    if (typeCount == idx2){
                        auto localOperator = tag_handler->get_op(positionPowers.at(nMaxVec[typeCount])[1].first);
                        localOperator *= positionPowers.at(nMaxVec[typeCount])[1].second;
                        local_op_vec.push_back(localOperator);
                    }
                    else {
                        local_op_vec.push_back(tag_handler->get_op(ident_vib.at(nMaxVec[idx2]).first));
                    }
                }
                typeCount++;
                ops.push_back(std::make_pair(local_op_vec, false));
                meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local,
                                                                      fillings_local, ops));
            }
            
        }
    }

    if(model.is_set("MEASURE[DisplacementSquared]")){
        int n_connectingmodes = model["vibronic_num_connectingmodes"].as<int>();
        int typeCount = 0; //typecount = type - 1 (for vibrational sites)
        for (std::size_t idx = 0; idx < n_particles_; idx++){
            for(std::size_t idx1 = 0; idx1 < n_vib_states_; idx1++){ 
                //if non existing lattice site: break. 
                //IMPORTANT: assumes underlying lattice sorting -> may be problematic
                //Assumed lattice sorting: intertwined with all connecting modes active. 
                //After an electronic site first come the local modes followed by the connecting modes
                if( (idx == n_particles_-1) && (idx1 >= (n_vib_states_-n_connectingmodes)) ) break; //VAL : break statement may be problematic... 
                std::string name = "DisplacementSquared"+std::to_string(idx)+"Mode"+std::to_string(idx1); 
                std::vector<pos_t> pos_internal(0);
                std::vector<std::vector<pos_t> > pos_local(0);
                pos_internal.push_back((idx1+1) + (idx*(n_vib_states_+n_ele_states_))); 
                maquis::cout << "pos of pushed back positions: " << idx1 << std::endl;
                pos_local.push_back(pos_internal);
                // Generates vector for the fillings and identity operators
                op_vec identities_local, fillings_local;
                for (std::size_t idx2 = 0; idx2 <= num_vibtypes; idx2++) {
                    identities_local.push_back(this->identity_matrix(idx2));
                    fillings_local.push_back(this->filling_matrix(idx2));
                }
                bond_element ops;
                op_vec local_op_vec;
                local_op_vec.push_back(tag_handler->get_op(ident_ele));
                for (std::size_t idx2 = 0; idx2 < num_vibtypes; idx2++) {
                    if (typeCount == idx2){
                        auto localOperator = tag_handler->get_op(positionPowers.at(nMaxVec[typeCount])[2].first);
                        localOperator *= positionPowers.at(nMaxVec[typeCount])[2].second;
                        local_op_vec.push_back(localOperator);
                    }
                    else{
                        local_op_vec.push_back(tag_handler->get_op(ident_vib.at(nMaxVec[idx2]).first));
                    }
                }
                typeCount++;
                ops.push_back(std::make_pair(local_op_vec, false));
                meas.push_back(new measurements::local_at<Matrix, U1>(name, lat, pos_local, identities_local,
                                                                      fillings_local, ops));
            }
            
        }
    }
    

        return meas;
    }


private:
    const Lattice& lat;
    value_type J_, epsilon_ ;
    BaseParameters& model;
    bool only_nn_;
    std::size_t L_, n_ele_states_, n_vib_states_, n_particles_;
    std::vector< Index<U1> > phys_indexes;
    std::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    std::unordered_map<int, std::pair<tag_type, value_type> > ident_vib;
    tag_type ident_ele, count_ele, count_ele_gs, create_ele, destroy_ele;
    /** Tag for the powers of the position/momentum operators */
    std::unordered_map<int, std::vector< std::pair<tag_type, value_type> >> positionPowers, momentumPowers;
    /** Maximum order of many-body coupling */
    int maxCoupling;
    std::vector<int> nMaxVec;
    int n_connectingmodes;
    int num_vibtypes; //number of different types for vibrational sites
    int check_maxVibMode;
};

#endif // DMRG_VIBRONIC