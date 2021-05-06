/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
 *               2020 by Robin Feldmann <robinfe@student.ethz.ch>
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

#ifndef MAQUIS_DMRG_NU1_TAG_GENERATION_HPP
#define MAQUIS_DMRG_NU1_TAG_GENERATION_HPP
#include <dmrg/models/lattice.h>
#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/matrix.hpp>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

#include <map>

using std::vector;
using std::pair;

template<class Matrix, class NU1>
class TagContainer {
    // Types definition
    typedef model_impl<Matrix, NU1> base;
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::op_t op_t;
    typedef typename std::vector<tag_type> operators_type;
    typedef typename base::measurements_type measurements_type;
    typedef typename Lattice::pos_t pos_t;
    typedef typename std::vector<pos_t> positions_type;
    typedef int term_type;
    typedef typename Matrix::value_type value_type;
    typedef typename NU1::charge charge_type;
    typedef typename std::pair<std::vector<std::size_t>, value_type> H_terms_type;

private:
    bool    verbose=false;
    std::size_t num_particle_types; // Number of particle types
    std::vector<unsigned int> vec_particles; // Number of particles per type
    std::vector<bool> isFermion; // Fermion 1, Boson 0
    // isFermion      = {1, 1, 0, 0, 1}
    std::vector<unsigned int> vec_orbitals; // Number of orbitals per type
    std::vector<unsigned int> vec_fer_bos; // vector that maps the particle types vector to a "fermion -- boson count vector"
    // vec_fer_bos    = {0, 1, 0, 1, 2}
    Lattice lat;
    std::vector<Index<NU1> > phys_indexes;
    std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler;
    std::vector<tag_type> fer_ident_tag, fer_filling_tag, fer_create_up_tag, fer_create_down_tag,
            fer_dest_up_tag, fer_dest_down_tag, bos_ident_tag,
            bos_create_tag, bos_dest_tag, bos_count_tag;


public:

    // Don't delete or the code doesn't compile... and don't ask me why
    TagContainer(){

    }

    TagContainer(const Lattice &lat_, std::shared_ptr<TagHandler<Matrix, NU1>> &tag_handler_)
            : lat(lat_), tag_handler(tag_handler_)
    {

        num_particle_types = boost::any_cast<std::size_t>(lat.get_parameter("num_particle_types"));
        vec_particles = boost::any_cast<std::vector<unsigned>>(lat.get_parameter("vec_particles"));
        isFermion = boost::any_cast<std::vector<bool>>(lat.get_parameter("isFermion"));
        vec_orbitals = boost::any_cast<std::vector<unsigned>>(lat.get_parameter("vec_orbitals"));
        vec_fer_bos = boost::any_cast<std::vector<unsigned>>(lat.get_parameter("vec_fer_bos"));

        //
        // Collect all fermions and all bosons in two different vectors.
        //
        std::size_t num_fer_types = 0; // Number of different fermions
        std::vector<std::size_t> fer_vec;
        std::size_t num_bos_types = 0; // Number of different bosons
        std::vector<std::size_t> bos_vec;
        for (std::size_t it = 0; it < num_particle_types; it++) {
            if (isFermion[it]) {
                num_fer_types++;
                fer_vec.push_back(vec_particles[it]);
            } else {
                num_bos_types++;
                bos_vec.push_back(vec_particles[it]);
            }
        }
        if (num_fer_types == 0 && num_bos_types == 0) throw std::runtime_error("No Fermions or Bosons?");

        std::vector<std::vector<charge_type> > fer_excited_states_array;
        std::vector<std::vector<charge_type> > bos_excited_states_array;
        charge_type empty_state(0);

        //
        //  filling the excited states arrays...
        //
        construct_excited_states(fer_excited_states_array, bos_excited_states_array);

        // All operators are collected in vectors
        // identity, creation operator (up, down), annihilation operator (up,down), occupation number operator.
        // Fermionic operators
        fer_ident_tag.reserve(num_fer_types);
        fer_filling_tag.reserve(num_fer_types);
        fer_create_up_tag.reserve(num_fer_types);
        fer_create_down_tag.reserve(num_fer_types);
        fer_dest_up_tag.reserve(num_fer_types);
        fer_dest_down_tag.reserve(num_fer_types);
        // Bosonic operators
        bos_ident_tag.reserve(num_bos_types);
        bos_create_tag.reserve(num_bos_types);
        bos_dest_tag.reserve(num_bos_types);
        bos_count_tag.reserve(num_bos_types);

        //
        // filling the tags...
        //
        create_tags(fer_excited_states_array, bos_excited_states_array, num_fer_types, num_bos_types);

        //
        // Next step is loading the types of physical indices.
        //
        phys_indexes.reserve(num_particle_types);

        if (num_fer_types > 0) {
            for (std::size_t it = 0; it < num_fer_types; it++) {
                Index<NU1> phys;
                phys.insert(std::make_pair(empty_state, 1));
                for (auto jt = fer_excited_states_array[it].begin(); jt < fer_excited_states_array[it].end(); jt++) {
                    phys.insert(std::make_pair(*jt, 1));
                }
                phys_indexes.push_back(phys);
            }
        }
        if (num_bos_types > 0) {
            for (std::size_t it = 0; it < num_bos_types; it++) {
                Index<NU1> phys;
                phys.insert(std::make_pair(empty_state, 1));
                for (auto jt = bos_excited_states_array[it].begin(); jt < bos_excited_states_array[it].end(); jt++) {
                    phys.insert(std::make_pair(*jt, 1));
                }
                phys_indexes.push_back(phys);
            }
        }


    }


    // +--------------------+
    // | REGISTER_ALL_TYPES |
    // +--------------------+
    /*! \brief
     *  This method is used to register all the operators contained in a given vector
     */
    std::vector<tag_type> register_all_types(std::vector<op_t> const &ops,   // -> vector of operators
                                             tag_detail::operator_kind kind)  // -> bosonic or fermionic
    {
        std::vector<tag_type> ret;
        for (std::size_t idx = 0; idx < ops.size(); idx++) {
            std::pair<tag_type, value_type> newtag = tag_handler->checked_register(ops[idx], kind);
            assert(newtag.first < tag_handler->size());
            assert(std::abs(newtag.second - value_type(1.)) == value_type());
            ret.push_back(newtag.first);
        }
        return ret;
    }

    // +--------------------------+
    // | CONSTRUCT EXCITED STATES |
    // +--------------------------+
    /*! \brief
     *  This function constructs the excited states array for fermions and bosons.
     * @param fer_excited_states_array
     * @param bos_excited_states_array
     */
    void construct_excited_states(std::vector<std::vector<charge_type>> &fer_excited_states_array,
                                  std::vector<std::vector<charge_type>> &bos_excited_states_array) {
        // Builds all the relevant charges.
        // The ones we are interested in are the empty state and each state in which
        // an orbital is populated
        charge_type empty_state(0);
        // Dimension of the excited_states array is the number of different particle types.
        // Now the excited states are populated; here it is important to distinguish between fermionic
        // and bosonic particles.
        // Fermions are associated with the basis <0,0>, <1,0>, <0,1>, and <1,1>,
        // whereas bosons are of type <0>, <1>, ..., <n>,<n+1>, where n is the number of bosons of the given type.
        // The unoccupied state is always implicit.
        for (std::size_t it = 0; it < num_particle_types; it++) {
            if (isFermion[it] == 1) {
                std::vector<charge_type> excited_states_it(3, 0);
                fer_excited_states_array.push_back(excited_states_it);
            } else if (isFermion[it] == 0) {
                std::vector<charge_type> excited_states_it(vec_particles[it] + 1, 0);
                bos_excited_states_array.push_back(excited_states_it);
            } else {
                throw std::runtime_error("Particle is not specified to be a boson or fermion.");
            }
        }
        // Prepare the excited states.
        // The fermion and boson counters are required for keeping track of the indices.
        int fcounter = 0;
        int bcounter = 0;
        for (std::size_t it = 0; it < num_particle_types; it++) {
            // Fermions: the same procedure for all fermions.
            if (isFermion[it] == 1) {
                fer_excited_states_array[fcounter][0][fcounter + it] = 1;
                fer_excited_states_array[fcounter][1][fcounter + it + 1] = 1;
                fer_excited_states_array[fcounter][2][fcounter + it] = 1;
                fer_excited_states_array[fcounter][2][fcounter + it + 1] = 1;
                fcounter++;
            }
                // Bosons: depends on the number of bosons of type "it"
            else if (isFermion[it] == 0) {
                for (std::size_t part_it = 0; part_it < vec_particles[it] + 1; part_it++) {
                    bos_excited_states_array[bcounter][part_it][it + fcounter] = part_it + 1;
                }
                bcounter++;
            }
        }
    } // construct_excited_states


    // +-------------+
    // | CREATE TAGS |
    // +-------------+
    /*! \brief
     * This function fills the tags for bosons and fermions.
     * The tags are needed to label the operators, which are later used for the construction of the MPO.
     * @param fer_excited_states_array
     * @param bos_excited_states_array
     * @param num_fer_types
     * @param num_bos_types
     */
    void create_tags(std::vector<std::vector<charge_type>> &fer_excited_states_array,
                     std::vector<std::vector<charge_type>> &bos_excited_states_array,
                     const std::size_t num_fer_types, const std::size_t num_bos_types) {

        std::vector<op_t> fer_ident, fer_filling, fer_create_up, fer_create_down, fer_dest_up, fer_dest_down;//, fer_count;
        std::vector<op_t> bos_ident, bos_create, bos_dest, bos_count;
        charge_type empty_state(0);
        //
        // loop over all fermion types to construct the opearators
        //
        int fcounter = 0;
        int bcounter = 0;
        for (std::size_t it = 0; it < num_particle_types; it++) {
            if (isFermion[it] == 1) {
                op_t ident, filling, create_up, create_down, destroy_up, destroy_down; //, count;
                // identities
                ident.insert_block(Matrix(1, 1, 1), empty_state, empty_state);
                for (auto jt = fer_excited_states_array[fcounter].begin();
                     jt < fer_excited_states_array[fcounter].end(); jt++)
                    ident.insert_block(Matrix(1, 1, 1), *jt, *jt);
                // filling
                filling.insert_block(Matrix(1, 1, 1), empty_state, empty_state);
                filling.insert_block(Matrix(1, 1, -1), fer_excited_states_array[fcounter][0],
                                     fer_excited_states_array[fcounter][0]);
                filling.insert_block(Matrix(1, 1, -1), fer_excited_states_array[fcounter][1],
                                     fer_excited_states_array[fcounter][1]);
                filling.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][2],
                                     fer_excited_states_array[fcounter][2]);
                // creators
                create_up.insert_block(Matrix(1, 1, 1), empty_state, fer_excited_states_array[fcounter][0]);
                create_up.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][1],
                                       fer_excited_states_array[fcounter][2]);
                create_down.insert_block(Matrix(1, 1, 1), empty_state, fer_excited_states_array[fcounter][1]);
                create_down.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][0],
                                         fer_excited_states_array[fcounter][2]);
                // destroyers
                destroy_up.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][0], empty_state);
                destroy_up.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][2],
                                        fer_excited_states_array[fcounter][1]);
                destroy_down.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][1], empty_state);
                destroy_down.insert_block(Matrix(1, 1, 1), fer_excited_states_array[fcounter][2],
                                          fer_excited_states_array[fcounter][0]);
                // push all the operators into the corresponding vectors.
                fer_ident.push_back(ident);
                fer_filling.push_back(filling);
                fer_create_up.push_back(create_up);
                fer_create_down.push_back(create_down);
                fer_dest_up.push_back(destroy_up);
                fer_dest_down.push_back(destroy_down);
                // fcounter
                fcounter++;
                if (this->verbose) {
                    std::cout << "Fermion ident" << std::endl;
                    std::cout << ident << std::endl;
                    std::cout << "Fermion filling" << std::endl;
                    std::cout << filling << std::endl;
                    std::cout << "Fermion create_up" << std::endl;
                    std::cout << create_up << std::endl;
                    std::cout << "Fermion create_down" << std::endl;
                    std::cout << create_down << std::endl;
                    std::cout << "Fermion destroy_up" << std::endl;
                    std::cout << destroy_up << std::endl;
                    std::cout << "Fermion destroy_down" << std::endl;
                    std::cout << destroy_down << std::endl;
                }
            } else if (isFermion[it] == 0) {
                op_t ident, create, destroy, count;
                // identities
                ident.insert_block(Matrix(1, 1, 1), empty_state, empty_state);
                for (auto jt = bos_excited_states_array[bcounter].begin();
                     jt < bos_excited_states_array[bcounter].end(); jt++)
                    ident.insert_block(Matrix(1, 1, 1), *jt, *jt);
                // creators & destroyers
                create.insert_block(Matrix(1, 1, 1), empty_state, bos_excited_states_array[bcounter][0]);
                destroy.insert_block(Matrix(1, 1, 1), bos_excited_states_array[bcounter][0], empty_state);
                for (std::size_t jt = 1; jt < vec_particles[it] + 1; jt++) {
                    // create
                    create.insert_block(Matrix(1, 1, std::sqrt(jt + 1)),
                                        bos_excited_states_array[bcounter][jt - 1],
                                        bos_excited_states_array[bcounter][jt]);
                    //destroy
                    destroy.insert_block(Matrix(1, 1, std::sqrt(jt + 1)),
                                         bos_excited_states_array[bcounter][jt],
                                         bos_excited_states_array[bcounter][jt - 1]);
                }
                // counters
                count.insert_block(Matrix(1, 1, 0), empty_state, empty_state);
                size_t jcont = 0;
                for (auto jt = bos_excited_states_array[bcounter].begin();
                     jt < bos_excited_states_array[bcounter].end(); jt++) {
                    jcont += 1;
                    count.insert_block(Matrix(1, 1, jcont), *jt, *jt);
                }
                // push all the operators into the corresponding vectors.
                bos_ident.push_back(ident);
                bos_count.push_back(count);
                bos_create.push_back(create);
                bos_dest.push_back(destroy);
                // fcounter
                bcounter++;
                if (this->verbose) {
                    std::cout << "Boson ident" << std::endl;
                    std::cout << ident << std::endl;
                    std::cout << "Boson count" << std::endl;
                    std::cout << count << std::endl;
                    std::cout << "Boson create" << std::endl;
                    std::cout << create << std::endl;
                    std::cout << "Boson destroy" << std::endl;
                    std::cout << destroy << std::endl;
                }
            }
        }
        //
        // Register all operators.
        //
        // Fermions
        if (num_fer_types > 0) {
            this->fer_ident_tag = register_all_types(fer_ident, tag_detail::bosonic);
            this->fer_filling_tag = register_all_types(fer_filling, tag_detail::bosonic);
            this->fer_create_up_tag = register_all_types(fer_create_up, tag_detail::fermionic);
            this->fer_create_down_tag = register_all_types(fer_create_down, tag_detail::fermionic);
            this->fer_dest_up_tag = register_all_types(fer_dest_up, tag_detail::fermionic);
            this->fer_dest_down_tag = register_all_types(fer_dest_down, tag_detail::fermionic);
        }
        // Bosons
        if (num_bos_types > 0) {
            this->bos_ident_tag = register_all_types(bos_ident, tag_detail::bosonic);
            this->bos_create_tag = register_all_types(bos_create, tag_detail::bosonic);
            this->bos_dest_tag = register_all_types(bos_dest, tag_detail::bosonic);
            this->bos_count_tag = register_all_types(bos_count, tag_detail::bosonic);
        }

    } //create_tags

    /*! \brief
     *  This functions assigns all variables which have been created upon construction of an instance of the class.
     */
    void
    get_all_variables(std::vector<Index<NU1>> &phys_indexes, std::shared_ptr<TagHandler<Matrix, NU1>> &tag_handler,
                      std::vector<tag_type> &fer_ident_tag, std::vector<tag_type> &fer_filling_tag,
                      std::vector<tag_type> &fer_create_up_tag, std::vector<tag_type> &fer_create_down_tag,
                      std::vector<tag_type> &fer_dest_up_tag, std::vector<tag_type> &fer_dest_down_tag,
                      std::vector<tag_type> &bos_ident_tag, std::vector<tag_type> &bos_create_tag,
                      std::vector<tag_type> &bos_dest_tag, std::vector<tag_type> &bos_count_tag)
    {
        phys_indexes        = this->phys_indexes;
        fer_ident_tag       = this->fer_ident_tag;
        fer_create_up_tag   = this->fer_create_up_tag;
        fer_dest_up_tag     = this->fer_dest_up_tag;
        bos_ident_tag       = this->bos_ident_tag;
        bos_dest_tag        = this->bos_dest_tag;
        fer_filling_tag     = this->fer_filling_tag;
        fer_create_down_tag = this->fer_create_down_tag;
        fer_dest_down_tag   = this->fer_dest_down_tag;
        bos_create_tag      = this->bos_create_tag;
        bos_count_tag       = this->bos_count_tag;
    } // get_all_variables
};
#endif //MAQUIS_DMRG_NU1_TAG_GENERATION_HPP
