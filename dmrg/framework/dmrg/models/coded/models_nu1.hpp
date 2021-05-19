/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

#ifndef MODELS_CODED_NU1_H
#define MODELS_CODED_NU1_H

/* internal includes */
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/prebo/prebo_TermGenerator.hpp"
#include "nu1_nBodyTerm.hpp"
#include "../prebo/prebo_parse_integrals.h"
#include "../model_helper.hpp"
#include "../measurements/prebo_particle_rdm.h"
/* external includes */
#include <algorithm>
#include <alps/hdf5/pair.hpp>
#include <boost/functional/hash.hpp>
#include <map>
#include <sstream>
#include <unordered_map>
#include <chrono>
// ========================
//  MODELS OF NU1 SYMMETRY
// ========================


// +----------------------------+
// | PRE BORN OPPENHEIMER MODEL |
// +----------------------------+
// Robin Feldmann
/**
 * @brief multicomponent DMRG calculation
 * @class preBO models_nu1.hpp
 * This class enables DMRG calculations with multiple indistinguishable fermionic (soon bosonic, too) particle types.
 * @tparam Matrix
 */
template<class Matrix>
class PreBO : public model_impl<Matrix, NU1> {
    // Types definition
    typedef model_impl<Matrix, NU1> base;
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    typedef typename base::term_descriptor term_descriptor;
    typedef typename std::vector<term_descriptor> terms_type;
    typedef typename base::op_t op_t;
    typedef typename std::vector<tag_type> operators_type;
    typedef typename base::measurements_type measurements_type;
    typedef typename Lattice::pos_t pos_t;                          // position on the lattice
    typedef typename Lattice::part_type part_type;                  // unsigned integer denoting the particle type
    typedef typename std::vector<pos_t> positions_type;
    typedef typename Matrix::value_type value_type;
    typedef typename NU1::charge charge_type;
    typedef typename std::pair<std::vector<std::size_t>, value_type> H_terms_type;

private:
    // +------------+
    // | Attributes |
    // +------------+
    const bool verbose = false;
    pos_t L;                        // Size of the DMRG lattice
    std::size_t num_particle_types; // Number of particle types
    std::vector<bool> isFermion;    // Fermion 1, Boson 0
    std::vector<int> vec_ini_state; // specify the initial state of the system.
    const Lattice& lat;
    BaseParameters& parms;
    std::vector<Index<NU1>> phys_indexes;
    std::shared_ptr<TagHandler<Matrix, NU1>> ptr_tag_handler;
    std::shared_ptr<prebo::TermGenerator<Matrix, NU1>> ptr_term_generator;
public:
    // +----------------+
    //  Main constructor
    // +----------------+
    PreBO(const Lattice& lat_, BaseParameters& parms_)
            : lat(lat_),
              parms(parms_),
              ptr_tag_handler(new table_type()),
              phys_indexes(0)
    {
        // Get only required variables from the lattice
        num_particle_types = lat.template get_prop<int>("num_particle_types");
        isFermion = lat.template get_prop<std::vector<bool>>("isFermion");
        vec_ini_state = lat.template get_prop<std::vector<int>>("vec_ini_state");
        L = lat.size();
        // Checks that the code has been compiled properly
        NU1* NU1_value = new NU1;
        int max_symm = NU1_value->get_dimension();
        delete NU1_value;
        // Checks whether the max symmetry number is sufficiently large
        unsigned int MIN = 0;
        unsigned int fcount = 0;
        for (auto it = isFermion.begin(); it < isFermion.end(); it++) {
            if (*it == 1)
                fcount++;
        }
        // The minimal number is twice the number of fermionic types plus the number
        // of bosonic types.
        MIN = 2 * fcount + num_particle_types - fcount;
        if (max_symm < MIN)
            throw std::runtime_error("Recompile QCMaquis with a larger value for DMRG_NUMSYMM");

        std::shared_ptr<Lattice> ptr_lat = std::make_shared<Lattice>(lat);
        ptr_term_generator = std::make_shared<prebo::TermGenerator<Matrix, NU1>>(ptr_lat, ptr_tag_handler);
        phys_indexes = ptr_term_generator->getPhysIndexes();

    } // Main constructor

    /**!
     * Create Hamiltonian terms
     * @return
     */
    void create_terms() {
        auto start = std::chrono::high_resolution_clock::now();
        std::pair<std::vector<chem::index_type<chem::Hamiltonian::PreBO>>, std::vector<double> > integrals = prebo::detail::parse_integrals<double, NU1>(parms, lat);
        this->terms_ = ptr_term_generator->generate_Hamiltonian(integrals);
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Construction of the Hamiltonian took "
                  << duration.count() << " seconds" << std::endl << std::endl;
    }


    //
    // +=========+
    // | METHODS |
    // +=========+
    //


    // +================+
    // | GETTER METHODS |
    // +================+
    /**
     * --> called in the constructor of the MPS Tensor Network
     * @param type
     * @return
     */
    Index<NU1> const& phys_dim(size_t type) const {
        return phys_indexes[type];
    }
    /**
     * --> Also called in the constructor of the MPS Tensor Network
     * The function recovers the identity or filling matrix.
     * @param type
     * @return
     */
    tag_type identity_matrix_tag(size_t type) const {
        auto vec_fer_bos = lat.template get_prop<std::vector<int>>("vec_fer_bos");
        // recover identity of the according particle type
        if (isFermion[type]) {
            auto tag_vec = ptr_term_generator->getTagContainer()->get_tag_vector(Type::Fermion, OpType::Ident, Spin::None);
            std::size_t fer = vec_fer_bos[type];
            return tag_vec[fer];
        }
        else {
            auto tag_vec = ptr_term_generator->getTagContainer()->get_tag_vector(Type::Boson, OpType::Ident, Spin::None);
            std::size_t bos = vec_fer_bos[type];
            return tag_vec[bos];
        }
    }
    /**
     * --> called in the constructor of the MPS Tensor Network
     * Returns the total quantum number associated to the Hamiltonian.
     * This means, the particle vector
     * @param parms
     * @return
     */
    typename NU1::charge total_quantum_numbers(BaseParameters& parms) const {
        // Not actually used.
        // int str = parms["PreBO_ParticleTypeVector"] ;
        // This is important:
        return NU1::charge(vec_ini_state);
    }
    /**
     * --> Also called in the constructor of the MPS Tensor Network
     * The function recovers the identity or filling matrix.
     * @param type
     * @return
     */
    tag_type filling_matrix_tag(size_t type) const {
        auto vec_fer_bos = lat.template get_prop<std::vector<int>>("vec_fer_bos");
        // recover identity of the according particle type
        if (isFermion[type]) {
            auto tag_vec = ptr_term_generator->getTagContainer()->get_tag_vector(Type::Fermion, OpType::Filling, Spin::None);
            std::size_t fer = vec_fer_bos[type];
            return tag_vec[fer];
        }
        else {
            auto tag_vec = ptr_term_generator->getTagContainer()->get_tag_vector(Type::Boson, OpType::Ident, Spin::None);
            std::size_t bos = vec_fer_bos[type];
            return tag_vec[bos];
        }
    }
    /**
     *
     * @param name
     * @param type
     * @return
     */
    tag_type get_operator_tag(std::string const& name, size_t type) const {
        throw std::runtime_error("Operator not valid for this model.");
    }

    void update(BaseParameters const& p) {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
    }
    //
    table_ptr operators_table() const {
        return ptr_tag_handler;
    }
    // TODO ALB Check if this really true, the index 0 has been hardcoded for the
    //         moment
    measurements_type measurements() const {
        measurements_type meas;

        std::vector<tag_type> identities;
        std::vector<tag_type> fillings;

        for (size_t p = 0; p <= lat.maximum_vertex_type(); ++p)
        {
            identities.push_back(this->identity_matrix_tag(p));
            fillings.push_back(this->filling_matrix_tag(p));
        }

        std::regex expression_1ParticleRDM("^MEASURE\\[1rdm\\]");
        std::smatch what;

        for (auto&& it: parms.get_range()) {
            std::string lhs = it.first;
            std::string name;
            // 1-RDM and transition-1RDM
            if (std::regex_match(lhs, what, expression_1ParticleRDM)) {
                name = "1ParticleRDM";
                meas.push_back( new measurements::PreBOParticleRDM<Matrix, NU1>(parms, lat, identities, fillings, ptr_term_generator));
            }
        }

        return meas;
    }
};

#endif
