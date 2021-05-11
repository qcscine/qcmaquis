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
#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"
#include "nu1_nBodyTerm.hpp"
#include "nu1_tag_generation.hpp"
//#include "integral_interface.h"
#include "../prebo/prebo_parse_integrals.h"
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
    pos_t L;                                 // Size of the DMRG lattice
    std::size_t num_particle_types;          // Number of particle types
    std::vector<unsigned int> vec_particles; // Number of particles per type
    std::vector<bool> isFermion;             // Fermion 1, Boson 0
    std::vector<unsigned int> vec_orbitals;  // number of orbitals for each particle type
    std::vector<unsigned int> vec_fer_bos;   // vector that maps the particle types vector. It counts fermions and bosons
    // to a "fermion -- boson count vector"
    // isFermion      = {1, 1, 0, 0, 1}
    // vec_fer_bos    = {0, 1, 0, 1, 2}
    std::vector<unsigned int> vec_ini_state; // specify the initial state of the system.
    const Lattice& lat;
    BaseParameters& model;
    std::vector<Index<NU1>> phys_indexes;
    std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler;
    std::vector<tag_type> fer_ident_tag, fer_filling_tag, fer_create_up_tag, fer_create_down_tag, fer_dest_up_tag,
            fer_dest_down_tag, bos_ident_tag, bos_create_tag, bos_dest_tag, bos_count_tag;
    terms_type terms_temp;
    // terms_type terms1RDM_;
    std::vector<pos_t> m_order;         // Ordering of the sites. By default starting from 0 to (L-1)
    std::vector<pos_t> m_inv_order;     // Inverse permutation of the order
public:
    // +----------------+
    //  Main constructor
    // +----------------+
    PreBO(const Lattice& lat_, BaseParameters& model_)
            : lat(lat_),
              model(model_),
              tag_handler(new table_type()),
              phys_indexes(0),
              num_particle_types(0),
              vec_particles(0),
              isFermion(0),
              vec_orbitals(0),
              vec_ini_state(0)
    // example H2O:
    // num_particle_types = 3
    // vec_particles = {10, 2, 1}
    // isFermion = {1, 1, 0}
    // vec_orbitals = {42, 42, 42}
    {
        //
        // Get all variables from the lattice
        //
        num_particle_types = boost::any_cast<std::size_t>(lat.get_prop_("num_particle_types"));
        vec_particles = boost::any_cast<std::vector<unsigned>>(lat.get_prop_("vec_particles"));
        isFermion = boost::any_cast<std::vector<bool>>(lat.get_prop_("isFermion"));
        vec_orbitals = boost::any_cast<std::vector<unsigned>>(lat.get_prop_("vec_orbitals"));
        vec_ini_state = boost::any_cast<std::vector<unsigned>>(lat.get_prop_("vec_ini_state"));
        vec_fer_bos = boost::any_cast<std::vector<unsigned>>(lat.get_prop_("vec_fer_bos"));
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

        TagContainer<Matrix, NU1> tag_container = TagContainer<Matrix, NU1>(lat, tag_handler);

        tag_container.get_all_variables(phys_indexes, tag_handler, fer_ident_tag, fer_filling_tag, fer_create_up_tag,
                                        fer_create_down_tag, fer_dest_up_tag, fer_dest_down_tag, bos_ident_tag,
                                        bos_create_tag, bos_dest_tag, bos_count_tag);
        // initialize the orbital lattice as the identity
        m_order.resize(L);
        m_inv_order.resize(L);
        for (pos_t i = 0; i < L; i++) {
            m_order[i] = i;
            m_inv_order[i] = i;
        }
        // If sites_order is given, permute the Hamiltonian MPO
        if (model.is_set("orbital_order")) {
            m_order = model["orbital_order"].as<std::vector<pos_t> >();
            if (m_order.size() != lat.size())
                throw std::runtime_error("orbital_order length is not the same as the number of orbitals\n");
            for (int p = 0; p < m_order.size(); ++p)
                m_inv_order[p] = std::distance(m_order.begin(), std::find(m_order.begin(), m_order.end(), p));
        }

    } // Main constructor

    /**!
     * Create Hamiltonian terms
     * @return
     */
    void create_terms() {
        auto start = std::chrono::high_resolution_clock::now();
        generate_Hamiltonian();
        auto stop = std::chrono::high_resolution_clock::now();
        auto duration = std::chrono::duration_cast<std::chrono::seconds>(stop - start);
        std::cout << "Construction of the Hamiltonian took "
                  << duration.count() << " seconds" << std::endl << std::endl;
    }

//    /**!
//     *
//     *
//     */
//    void generate_terms1RDM(int i, int mu, int nu) override {
//        //
//        // < Psi | a+^i_mu a^j_nu | Psi> => Generate operator string for particle type i and sites mu and nu
//        //
//        this->terms1RDM_.resize(0);
//        // the first element corresponds to the particle type, the second one is
//        // the orbital index.
//        std::vector<std::pair<part_type, pos_t>> nbody_term;
//        nbody_term.push_back(std::make_pair(i, mu));
//        nbody_term.push_back(std::make_pair(i, nu));
//
//        // Assert that the particle type index is less than the number of particle
//        // types and that the orbital index is less than the number of orbitals of
//        // the given type
//        for (auto const& iter : nbody_term) {
//            assert(iter.first < num_particle_types);
//            assert(iter.second < vec_orbitals.at(iter.first));
//        }
//
//        // The nbody term is now ready to be transformed to the second quantized
//        // operators with proper symmetry and spin configurations.
//        NBodyTerm nBodyTerm = NBodyTerm(nbody_term, isFermion, vec_orbitals, m_inv_order);
//
//        std::vector<std::vector<SymbolicOperator>> vec_SymOpStr = nBodyTerm.getVecSymOpStr();
//
//        for (auto const& OpStr : vec_SymOpStr) {
//            positions_type pos;
//            operators_type ops;
//            Symbols2Tag(pos, ops, OpStr);
//            value_type scaling = 1.;
//            std::pair<term_descriptor, bool> ret = arrange_operators(pos, ops, scaling, tag_handler);
//            if (!ret.second) {
//                auto term = ret.first;
//                term.coeff = scaling;
//                this->terms1RDM_.push_back(term);
//                if (this->verbose) {
//                    std::cout << term << std::endl;
//                }
//            }
//        }
//    }

    //
    // +=========+
    // | METHODS |
    // +=========+
    //
    // +-------------------+
    // | CLEAN_HAMILTONIAN |
    // +-------------------+
    /**! \brief
     * The purpose of this function is to remove all duplicates of terms that
     * appear more than one time in the Hamiltonian and add the coefficient if
     * they belong to the same combination of operators.
     *     => This drastically reduces the bond dimensions.
     */
    void clean_hamiltonian() {
        std::cout << "======================================================" << std::endl;
        std::cout << "Cleaning the Hamiltonian..." << std::endl;
        std::unordered_map<term_descriptor, value_type, term_descriptor_hasher> count_map;

        for (auto const& term1 : this->terms_temp) {
            // If key not found in map iterator to end is returned
            if (count_map.find(term1) == count_map.end()) {
                count_map[term1] = term1.coeff;
            }
                // If key found then iterator to that key is returned
            else {
                count_map[term1] += term1.coeff;
            }
        }

        std::cout << "======================================================" << std::endl;
        std::cout << "The final size of the Hamiltonian is:" << std::endl << std::endl;
        std::cout << count_map.size() << std::endl << std::endl;
        std::cout << "======================================================" << std::endl;

        this->terms_.reserve(count_map.size());
        for (const auto& iter : count_map) {
            auto term = iter.first;
            term.coeff = iter.second;
            this->terms_.push_back(std::move(term));
        }

    } // clean_hamiltonian

    struct term_descriptor_hasher {
        std::size_t operator()(const term_descriptor& key) const {
            using boost::hash_combine;
            using boost::hash_value;

            // Start with a hash value of 0    .
            std::size_t seed = 0;

            // Modify 'seed' by XORing and bit-shifting in
            // one member of 'Key' after the other:
            for (unsigned int i = 0; i < key.size(); i++) {
                hash_combine(seed, key.position(i));
                hash_combine(seed, key.operator_tag(i));
            }
            // Return the result.
            return seed;
        }
    };

    /**
     * This method takes the symbolic operator string and populates operators and
     * positions vector.
     * @param pos
     * @param ops
     * @param SymOpStr
     */
    void Symbols2Tag(positions_type& pos, operators_type& ops, const std::vector<SymbolicOperator>& SymOpStr) {
        for (auto const& SymOp : SymOpStr) {
            pos.push_back(SymOp.getSite());
            if (SymOp.getOpType() == SymbolicOperator::FILLING) {
                ops.push_back(fer_filling_tag[vec_fer_bos[SymOp.getPartType()]]);
            }
            else if (SymOp.getOpType() == SymbolicOperator::IDENT) {
                if (isFermion[SymOp.getPartType()])
                    ops.push_back(fer_ident_tag[vec_fer_bos[SymOp.getPartType()]]);
                else
                    ops.push_back(bos_ident_tag[vec_fer_bos[SymOp.getPartType()]]);
            }
            else if (SymOp.getOpType() == SymbolicOperator::CREATE) {
                if (SymOp.getSpin() == SymbolicOperator::DOWN)
                    ops.push_back(fer_create_down_tag[vec_fer_bos[SymOp.getPartType()]]);
                else if (SymOp.getSpin() == SymbolicOperator::UP)
                    ops.push_back(fer_create_up_tag[vec_fer_bos[SymOp.getPartType()]]);
                else if (SymOp.getSpin() == SymbolicOperator::ZERO)
                    ops.push_back(bos_create_tag[vec_fer_bos[SymOp.getPartType()]]);
            }
            else if (SymOp.getOpType() == SymbolicOperator::ANNIHILATE) {
                if (SymOp.getSpin() == SymbolicOperator::DOWN)
                    ops.push_back(fer_dest_down_tag[vec_fer_bos[SymOp.getPartType()]]);
                else if (SymOp.getSpin() == SymbolicOperator::UP)
                    ops.push_back(fer_dest_up_tag[vec_fer_bos[SymOp.getPartType()]]);
                else if (SymOp.getSpin() == SymbolicOperator::ZERO)
                    ops.push_back(bos_dest_tag[vec_fer_bos[SymOp.getPartType()]]);
            }
        }
    }

    // +-------------------+
    // | ARRANGE_OPERATORS |
    // +-------------------+
    /**
     * This routine is used to take a list of SQ operators and get the
     * corresponding list of tags. Operators centered on the same center are
     * merged together.
     * @param positions
     * @param operators
     * @param tag_handler
     * @return
     */
    std::pair<term_descriptor, bool> arrange_operators(positions_type const& positions, operators_type const& operators,
                                                       value_type& scaling,
                                                       std::shared_ptr<TagHandler<Matrix, NU1>> tag_handler) {
        // Safety check
        assert(positions.size() == operators.size());
        bool FoundZero = false;
        // Types definition
        typedef Lattice::pos_t pos_t;
        typedef typename Matrix::value_type value_type;
        typedef typename OPTable<Matrix, NU1>::tag_type tag_type;
        typedef std::pair<pos_t, tag_type> pos_op_t;
        // Variables definition
        term_descriptor term;
        std::vector<pos_op_t> pos_ops;
        std::transform(positions.begin(), positions.end(), operators.begin(), std::back_inserter(pos_ops),
                       std::make_pair<pos_t const&, tag_type const&>);
        std::stable_sort(pos_ops.begin(), pos_ops.end(), generate_mpo::compare<pos_op_t>);
        // Now that the operators are properly sorted, the formation of the tag can
        // start. Note that get_product_tag returns a new tag if the product does
        // not exists, otherwise returns the existing tag
        for (size_t opnr = 0; opnr < pos_ops.size();) {
            tag_type product = pos_ops[opnr].second;
            size_t range_end = opnr + 1;
            while (range_end < pos_ops.size() && pos_ops[range_end].first == pos_ops[opnr].first) {
                value_type scale = 1.0;
                if (tag_handler->is_null(pos_ops[range_end].second, product))
                    FoundZero = true;
                boost::tie(product, scale) = tag_handler->get_product_tag(pos_ops[range_end].second, product);
                scaling *= scale;
                range_end++;
            }
            term.push_back(boost::make_tuple(pos_ops[opnr].first, product));
            opnr = range_end;
        }
        // std::cout << "Overall scaling" << std::endl;
        // std::cout << scaling << std::endl;
        return std::make_pair(term, FoundZero);
    } // arrange_operators

    // +----------+
    // | ADD TERM |
    // +----------+
    /**
     * Adds a single term to the Hamiltonian object
     * @param positions
     * @param operators
     * @param coeff
     */
    void add_term(positions_type const& positions, operators_type const& operators, value_type const& coeff) {
        static int count = 0;
        value_type scaling = 1.;
        std::pair<term_descriptor, bool> ret = arrange_operators(positions, operators, scaling, tag_handler);
        if (!ret.second) {
            count++;
            auto term = ret.first;
            term.coeff = coeff * scaling;
            this->terms_temp.push_back(term);
            if (this->verbose) {
                std::cout << term << std::endl;
                std::cout << "Operator count = " << count << std::endl;
            }
        }
    } //

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
        // recover identity of the according particle type
        if (isFermion[type]) {
            std::size_t fer = vec_fer_bos[type];
            return fer_ident_tag[fer];
        }
        else {
            std::size_t bos = vec_fer_bos[type];
            return bos_ident_tag[bos];
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
        // recover filling of the according particle type
        if (isFermion[type]) {
            std::size_t fer = vec_fer_bos[type];
            return fer_filling_tag[fer];
        }
        else {
            std::size_t bos = vec_fer_bos[type];
            return bos_ident_tag[bos];
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
        return tag_handler;
    }
    // TODO ALB Check if this really true, the index 0 has been hardcoded for the
    //         moment
    measurements_type measurements() const {
        // typedef std::vector<op_t> op_vec;
        // typedef std::vector<std::pair<op_vec, bool> > bond_element;
        measurements_type meas;
        // if (model["MEASURE[Density]"]) {
        //    std::string name = "Density";
        //    meas.push_back( new measurements::average<Matrix, NU1>(name, lat,
        //                                                           tag_handler->get_ops(ident),
        //                                                           tag_handler->get_ops(ident),
        //                                                           tag_handler->get_ops(count))
        //                                                           ) ;
        //}
        // if (model["MEASURE[Local density]"]) {
        //    std::string name = "Local density";
        //    meas.push_back( new measurements::local<Matrix, NU1>(name, lat,
        //                                                         tag_handler->get_ops(ident),
        //                                                         tag_handler->get_ops(ident),
        //                                                         tag_handler->get_ops(count))
        //                                                         ) ;
        //}
        // if (model["MEASURE[Onebody density matrix]"]) {
        //    std::string name = "Onebody density matrix";
        //    bond_element ops;
        //    ops.push_back( std::make_pair( tag_handler->get_ops(create), false) );
        //    ops.push_back( std::make_pair( tag_handler->get_ops(destroy), false)
        //    ); meas.push_back( new measurements::correlations<Matrix, NU1>(name,
        //    lat,
        //                                                   op_vec(1,this->identity_matrix(0)),
        //                                                   op_vec(1,this->filling_matrix(0)),
        //											       ops, true, false) )
        //;
        //}
        return meas;
    }

    void generate_Hamiltonian() {
        // +-------------------------------------+
        // | New Construction of the Hamiltonian |
        // +-------------------------------------+
        std::cout << "======================================================" << std::endl;
        std::cout << "Starting to construct the Hamiltonian..." << std::endl;
        // ORDER BEGIN
        bool sorting = false;
        // Load ordering and determine inverse ordering
        // Open the integral file and do a loop over it.
        std::string integral_file = model["integral_file"];
        if (!boost::filesystem::exists(integral_file))
            throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
        std::ifstream orb_file;
        std::string line_string;
        //
        orb_file.open(integral_file.c_str());

        // +-------------------------------------+
        // | New Construction of the Hamiltonian |
        // +-------------------------------------+
        std::vector<std::pair<part_type, pos_t>> nbody_term;

        std::pair<std::vector<chem::index_type<chem::Hamiltonian::PreBO>>, std::vector<double> > integrals = prebo::detail::parse_integrals<double, NU1>(model, lat);

        auto getTermOrder = [] (const chem::index_type<chem::Hamiltonian::PreBO>& indices) {
            auto ans = std::count(indices.begin(), indices.end(), -1);
            if (ans == 8)
                return 0;
            else if (ans == 4)
                return 1;
            else if (ans == 8)
                return 2;
            else {
                throw std::runtime_error("Term has invlaid order");
            }
        };

        unsigned termOrder=0;
        for (auto line=0; line<integrals.first.size(); ++line) {
            /// -- Initialization --
            nbody_term.clear();
            // Extract the order of the term:
            termOrder = getTermOrder(integrals.first[line]);
            // Nuclear repulsion:
            const auto& indices = integrals.first[line];
            const auto& integral = integrals.second[line];
            if (termOrder==0) {
                positions_type pos{0};
                operators_type ops;
                if (isFermion[lat.template get_prop<int>("type", pos)])
                    ops.push_back(fer_ident_tag[lat.template get_prop<int>("type", pos)]);
                else
                    ops.push_back(bos_ident_tag[lat.template get_prop<int>("type", pos)]);
                add_term(pos, ops, integral);
                continue;
            }
            // 1-body and 2-body terms
            for (auto i=0; i<termOrder*4; i+=2) {
                nbody_term.push_back(std::make_pair(indices[i], indices[i+1]));
            }
            // Assert tat for a one-body term, the particle types are equal
            // and that for a two-body term, the rule a(i)+ a(m)+ a(m) a(i) is
            // followed.
            assert((termOrder == 1 && nbody_term[0].first == nbody_term[1].first) ||
                   (termOrder == 2 && nbody_term[0].first == nbody_term[3].first &&
                    nbody_term[1].first == nbody_term[2].first));
            // Assert that the particle type index is less than the number of particle
            // types and that the orbital index is less than the number of orbitals of
            // the given type
            #ifndef NDEBUG
            for (auto const& iter : nbody_term) {
                assert(iter.first < num_particle_types);
                assert(iter.second < vec_orbitals.at(iter.first));
            }
            #endif
            // The nbody term is now ready to be transformed to the second quantized
            // operators with proper symmetry and spin configurations.
            NBodyTerm nBodyTerm = NBodyTerm(nbody_term, isFermion, vec_orbitals, m_inv_order);

            std::vector<std::vector<SymbolicOperator>> vec_SymOpStr = nBodyTerm.getVecSymOpStr();

            for (auto const& OpStr : vec_SymOpStr) {
                positions_type pos;
                operators_type ops;
                Symbols2Tag(pos, ops, OpStr);
                add_term(pos, ops, integral);
            }

        }
        //while (getline(orb_file, line_string)) {
        //    if (line_string[0] == '#' || line_string=="")
        //        continue;
        //    // initialize integral value
        //    value_type integral;
        //    // initialize splitted line
        //    std::vector<std::string> line_splitted;
        //    // -- Main data parsing --
        //    // Trim leading and final spaces in the string.
        //    line_string.erase(line_string.begin(),
        //                      std::find_if(line_string.begin(), line_string.end(), [&](int ch) { return !std::isspace(ch); }));
        //    line_string.erase(
        //            std::find_if(line_string.rbegin(), line_string.rend(), [&](int ch) { return !std::isspace(ch); }).base(),
        //            line_string.end());
        //    // Split the string
        //    boost::split(line_splitted, line_string, boost::is_any_of(" "), boost::token_compress_on);
        //    // EMPTY LINE
        //    // if (line_splitted.size() == 0) continue;

        //    // Last value in string is assigned to the integral value
        //    integral = atof(line_splitted[line_splitted.size() - 1].c_str());

        //    //
        //    // Nuclear repulsion
        //    //
        //    if (line_splitted.size() == 1 ) {
        //        positions_type pos{0};
        //        operators_type ops;
        //        if (isFermion[lat.template get_prop<int>("type", pos)])
        //            ops.push_back(fer_ident_tag[lat.template get_prop<int>("type", pos)]);
        //        else
        //            ops.push_back(bos_ident_tag[lat.template get_prop<int>("type", pos)]);
        //        add_term(pos, ops, integral);
        //        continue;
        //    }
        //    // remove integral value from vector
        //    // now the vector contains all 2nd quantization operators.
        //    line_splitted.pop_back();
        //    assert(line_splitted.size() == 2 || line_splitted.size() == 4);
        //    // loop over all 2nd quant. operators in vector.
        //    for (const auto& sq_op_str : line_splitted) {
        //        std::vector<std::string> temp;
        //        boost::split(temp, sq_op_str, boost::is_any_of("-"));
        //        assert(temp.size() == 2);
        //        // the first element corresponds to the particle type, the second one is
        //        // the orbital index.
        //        nbody_term.push_back(std::make_pair(std::stoul(temp[0]), std::stoul(temp[1])));
        //    }
        //    // Assert tat for a one-body term, the particle types are equal
        //    // and that for a two-body term, the rule a(i)+ a(m)+ a(m) a(i) is
        //    // followed.
        //    assert((line_splitted.size() == 2 && nbody_term[0].first == nbody_term[1].first) ||
        //           (line_splitted.size() == 4 && nbody_term[0].first == nbody_term[3].first &&
        //            nbody_term[1].first == nbody_term[2].first));
        //    // Assert that the particle type index is less than the number of particle
        //    // types and that the orbital index is less than the number of orbitals of
        //    // the given type
        //    for (auto const& iter : nbody_term) {
        //        assert(iter.first < num_particle_types);
        //        assert(iter.second < vec_orbitals.at(iter.first));
        //    }

        //    // The nbody term is now ready to be transformed to the second quantized
        //    // operators with proper symmetry and spin configurations.
        //    NBodyTerm nBodyTerm = NBodyTerm(nbody_term, isFermion, vec_orbitals, m_inv_order);

        //    std::vector<std::vector<SymbolicOperator>> vec_SymOpStr = nBodyTerm.getVecSymOpStr();

        //    for (auto const& OpStr : vec_SymOpStr) {
        //        positions_type pos;
        //        operators_type ops;
        //        Symbols2Tag(pos, ops, OpStr);
        //        add_term(pos, ops, integral);
        //    }
        //}

        //
        // And finally ... tidy up the mess.
        //
        clean_hamiltonian();
    }

    std::vector<std::pair<value_type, std::vector<SymbolicOperator>>> multiply(
            const std::vector<std::pair<value_type, std::vector<SymbolicOperator>>> & lhs,
            const std::vector<std::pair<value_type, std::vector<SymbolicOperator>>> & rhs)
    {
        std::vector<std::pair<value_type, std::vector<SymbolicOperator>>> res;
        res.reserve(lhs.size()*rhs.size());
        // Multiply the two strings of symbolic operators.
        // Loop over terms on the left hand side and multiply for each term with every other term on the right hand side
        for (auto const & itlhs: lhs)
            for (auto const & itrhs: rhs) {
                value_type val=itlhs.first*itrhs.first;
                std::vector<SymbolicOperator> tmp_str = itlhs.second;
                // Concatenate two strings of operators
                tmp_str.insert(tmp_str.end(), itrhs.second.begin(), itrhs.second.end());
                res.push_back(std::make_pair(val, tmp_str));
            }
        return res;
    }

//    std::vector<std::pair<value_type, std::vector<SymbolicOperator>>>
//    generate_transitionOps(const part_type i, const pos_t mu, const int which) {
//
//        std::vector<std::pair<value_type, std::vector<SymbolicOperator>>> transition_op;
//        std::vector<SymbolicOperator> SymOp_temp;
//        positions_type pos{mu};
//        operators_type ops;
//        value_type scaling = 1.;
//
//        if (isFermion[i]) {
//            switch (which) {
//                case 1:
//                    //
//                    // Operator 1: 1 - n_up - n_down + n_up n_down
//                    //
//                    // 1
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::IDENT, i));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // + n_up n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 2:
//                    //
//                    // Operator 2: c_down  - n_up c_down
//                    //
//                    // F c_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_up c_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 3:
//                    //
//                    // Operator 3: c_up - n_down c_up
//                    //
//                    // c_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_down c_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 4:
//                    //
//                    // Operator 4:  c_down c_up
//                    //
//                    // c_down c_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 5:
//                    //
//                    // Operator 5: c+_down  - n_up c+_down
//                    //
//                    // c+_down F
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_up c+_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 6:
//                    //
//                    // Operator 6: n_down - n_up n_down
//                    //
//                    // n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_up n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 7:
//                    //
//                    // Operator 7:  c+_down c_up
//                    //
//                    transition_op.resize(0);
//                    // c+_down F c_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 8:
//                    //
//                    // Operator 8: - n_down c_up
//                    //
//                    transition_op.resize(0);
//                    // - n_down c_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 9:
//                    //
//                    // Operator 9: c+_up - n_down c+_up
//                    //
//                    // c+_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_down c+_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 10:
//                    //
//                    // Operator 10:  c_down c+_up
//                    //
//                    // c_down c+_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 11:
//                    //
//                    // Operator 11: n_up - n_up n_down
//                    //
//                    // n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    // - n_up n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 12:
//                    //
//                    // Operator 12: n_up c_down
//                    //
//                    // n_up c_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 13:
//                    //
//                    // Operator 13:  c+_down c+_up
//                    //
//                    // c+_down c+_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 14:
//                    //
//                    // Operator 14: - n_down c+_up
//                    //
//                    // - n_down c+_up
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    scaling = -1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 15:
//                    //
//                    // Operator 15: n_up c+_down
//                    //
//                    // n_up c+_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//                case 16:
//                    //
//                    // Operator 16: n_up n_down
//                    //
//                    // n_up n_down
//                    SymOp_temp.clear();
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::UP));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::CREATE, i, SymbolicOperator::DOWN));
//                    SymOp_temp.push_back(SymbolicOperator(mu, SymbolicOperator::ANNIHILATE, i, SymbolicOperator::DOWN));
//                    scaling = 1.;
//                    transition_op.push_back(std::make_pair(scaling, SymOp_temp));
//                    break;
//            }
//        }
//        else {
//            std::cout << "Bosonic Transition Operators not implemented, yet." << std::endl;
//        }
//        return transition_op;
//    }
//
//
//    std::vector<std::pair<value_type, std::vector<SymbolicOperator>>>
//    generate_transitionOps(const part_type i, const part_type j, const pos_t mu, const pos_t nu, const int which1, const int which2) {
//
//        SymbolicJordanWigner JW;
//
//        std::vector<std::pair<value_type,std::vector<SymbolicOperator>>> transition_ops1 =
//                generate_transitionOps(i, mu, which1);
//        std::vector<std::pair<value_type,std::vector<SymbolicOperator>>> transition_ops2 =
//                generate_transitionOps(j, nu, which2);
//
//        if (i!=j) {
//            for (size_t i=0; i<transition_ops1.size(); i++) {
//                JW = SymbolicJordanWigner(transition_ops1.at(i).second);
//                transition_ops1.at(i).second = JW.getSymOpStr();
//            }
//            for (size_t i=0; i<transition_ops2.size(); i++) {
//                JW = SymbolicJordanWigner(transition_ops2.at(i).second);
//                transition_ops2.at(i).second = JW.getSymOpStr();
//            }
//        }
//        std::vector<std::pair<value_type,std::vector<SymbolicOperator>>> transition_ops_res = multiply(transition_ops1,
//                                                                                                       transition_ops2);
//        if (i==j) {
//            for (size_t i=0; i<transition_ops_res.size(); i++) {
//                JW = SymbolicJordanWigner(transition_ops_res.at(i).second);
//                transition_ops_res.at(i).second = JW.getSymOpStr();
//            }
//        }
//        return transition_ops_res;
//
//    }
//
//    void generate_termsTransitionOp(int i, int mu, int which) {
//
//        this->terms1oRDM_.resize(0);
//
//        std::vector<std::pair<value_type,std::vector<SymbolicOperator>>> transition_ops =
//                generate_transitionOps(i, mu, which);
//
//        SymbolicJordanWigner JW;
//        for (size_t i=0; i<transition_ops.size(); i++) {
//            JW = SymbolicJordanWigner(transition_ops.at(i).second);
//            transition_ops.at(i).second = JW.getSymOpStr();
//        }
//
//        for (auto const& transOp : transition_ops ) {
//            positions_type pos;
//            operators_type ops;
//            Symbols2Tag(pos, ops, transOp.second);
//            value_type scaling = transOp.first;
//            std::pair<term_descriptor, bool> ret = arrange_operators(pos, ops, scaling, tag_handler);
//            if (!ret.second) {
//                auto term = ret.first;
//                term.coeff = scaling;
//                this->terms1oRDM_.push_back(term);
//                if (this->verbose) {
//                    std::cout << term << std::endl;
//                }
//            }
//        }
//
//    }
//
//
//    void generate_termsTransitionOp(int i, int j, int mu, int nu, int which1, int which2) {
//
//        this->terms2oRDM_.resize(0);
//
//        std::vector<std::pair<value_type,std::vector<SymbolicOperator>>> transition_ops =
//                generate_transitionOps(i, j, mu, nu, which1, which2);
//
//        for (auto const& transOp : transition_ops ) {
//            positions_type pos;
//            operators_type ops;
//            Symbols2Tag(pos, ops, transOp.second);
//            value_type scaling = transOp.first;
//            std::pair<term_descriptor, bool> ret = arrange_operators(pos, ops, scaling, tag_handler);
//            if (!ret.second) {
//                auto term = ret.first;
//                term.coeff = scaling;
//                this->terms2oRDM_.push_back(term);
//                if (this->verbose) {
//                    std::cout << term << std::endl;
//                }
//            }
//        }
//    }

};

#endif
