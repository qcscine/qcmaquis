/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef GENERATE_MPO_TAGGED_MPO_MAKER_H
#define GENERATE_MPO_TAGGED_MPO_MAKER_H

#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/block_matrix/block_matrix.h"
#include "dmrg/block_matrix/block_matrix_algorithms.h"
#include "dmrg/block_matrix/symmetry.h"
#include "dmrg/mp_tensors/mpo.h"
#include "dmrg/mp_tensors/mpo_ops.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/model.h"
#include <string>
#include <sstream>
#include <tuple>

namespace generate_mpo {

namespace detail {

/**
 * @brief Structure representing the key of a prempo object
 *
 * This class represents the objects that we assign to each "b" value
 * of the MPO.
 */
template <typename pos_t, typename tag_type, typename index_type>
struct prempo_key {
    typedef std::pair<pos_t, tag_type> pos_op_type;
    enum kind_type {trivial_left, bulk, bulk_no_merge, trivial_right};
    // Class members
    kind_type kind;
    std::vector<pos_op_type> pos_op;
    index_type offset;
    /** @brief Default constructor */
    prempo_key(kind_type k_=bulk, index_type o_=0) : kind(k_), offset(o_) { }

    /** @brief Constructor from a vector of position/tag pairs */
    prempo_key(std::vector<pos_op_type> const& po_, index_type o_=0) : kind(bulk), pos_op(po_), offset(o_) { }

    bool operator==(prempo_key const& lhs) const {
        if (kind != lhs.kind)
            return false;
        if (kind == trivial_left)
            return true;
        if (kind == trivial_right)
            return true;
        return (pos_op == lhs.pos_op) && (offset == lhs.offset);
    }

    bool operator<(prempo_key const& lhs) const {
        if (kind != lhs.kind) return kind < lhs.kind;
        //if (pos_op.size() != lhs.pos_op.size()) return pos_op.size() < lhs.pos_op.size();
        return (pos_op == lhs.pos_op) ? offset < lhs.offset : pos_op < lhs.pos_op;
    }
};

} // namespace detail

template <typename pos_t, typename tag_type, typename index_type>
std::ostream& operator << (std::ostream& os, detail::prempo_key<pos_t, tag_type, index_type> key) {
    auto s = key.pos_op.size();
    for (int i = 0; i < s; ++i)
        os << key.pos_op[i].first << ":" << key.pos_op[i].second << ", ";
    os << "o" << key.offset;
    return os;
}

/**
 * @brief TaggedMPOMaker class
 *
 * Intermediate class that takes the terms from a model and generates an object that can be
 * then fed to the MPO constructor
 * The class has the following attributes:
 *
 * lat: reference lattice
 * identities/identities_full: tag associated with the identity operator
 * fillings: tag associated with the filling operator
 * lenght: lattice size
 * tag_handler: map in which all the operators are stored as tags
 * prempo: vector of prempo_map_type objects, which are map associating pairs of prempo_key objects to values
 *         (i.e., coefficients of the Hamiltonian). The pair of objects, which are the keyword for the dictionary,
 *         are the labels which are associated to the MPO bond, i.e. the operators which have been
 *         applied before and/or after the current site
 * trivial_left/trivial_right: operators to be put on the left/right of each strings of SQ operators
 * site_terms: map storing operators acting ONLY on the site (diagonal one-body operators)
 * leftmost_right/rightmost_left: last and first site where there are no trivial operators
 * core_energy: coefficient of the full identity operator
 */
template<class Matrix, class SymmGroup>
class TaggedMPOMaker
{
    using scale_type = typename Matrix::value_type;
    using index_type = typename MPOTensor<Matrix, SymmGroup>::index_type;
    using op_t = typename OPTable<Matrix, SymmGroup>::op_t;
    using pos_t = Lattice::pos_t;
    using tag_type = typename OperatorTagTerm<Matrix, SymmGroup>::tag_type;
    using pos_op_type = typename OperatorTagTerm<Matrix, SymmGroup>::op_pair_t;
    using tag_block = boost::tuple<std::size_t, std::size_t, tag_type, scale_type>;
    using term_descriptor = ::term_descriptor<typename Matrix::value_type>;
    using tag_vec = std::vector<tag_type>;
    using prempo_key_type = detail::prempo_key<pos_t, tag_type, index_type>;
    using prempo_value_type = std::pair<tag_type, scale_type>;
    // TODO: consider moving to hashmap
    using prempo_map_type =  std::multimap<std::pair<prempo_key_type, prempo_key_type>, prempo_value_type,
                                           compare_pair_inverse<std::pair<prempo_key_type, prempo_key_type> > >;
    enum merge_kind {attach, detach};

public:
    /** @brief Class constructor from a lattice and a model */
    TaggedMPOMaker(Lattice const& lat_, Model<Matrix,SymmGroup> const& model)
        : lat(lat_), length(lat.size()), tag_handler(model.operators_table()), prempo(length),
          trivial_left(prempo_key_type::trivial_left), trivial_right(prempo_key_type::trivial_right),
          leftmost_right(length), rightmost_left(0), finalized(false), verbose(true), core_energy(0.)
    {
        for (int iType = 0; iType < lat.getMaxType(); iType++) {
            identities.push_back(model.identity_matrix_tag(iType));
            fillings.push_back(model.filling_matrix_tag(iType));
            try {
                identities_full.push_back(model.get_operator_tag("ident_full", iType));
            }
            catch (std::runtime_error const & e) {}
        }
        for (const auto& iTerm: model.hamiltonian_terms())
            this->add_term(iTerm);
    }

    /**
     * @brief Class constructor from a lattice and a vector of terms.
     * This method should be used for constructing operators != from the Hamiltonian.
     */
    TaggedMPOMaker(Lattice const& lat_, tag_vec const & i_, tag_vec const & i_f_, tag_vec const & f_,
                   std::shared_ptr<TagHandler<Matrix, SymmGroup> > th_, typename Model<Matrix, SymmGroup>::terms_type const& terms)
    : lat(lat_), identities(i_), identities_full(i_f_), fillings(f_), length(lat.size()), tag_handler(th_), prempo(length),
      trivial_left(prempo_key_type::trivial_left), trivial_right(prempo_key_type::trivial_right), leftmost_right(length),
      rightmost_left(0), finalized(false), verbose(false), core_energy(0.)
    {
        //for (size_t p = 0; p < length-1; ++p)
        //    prempo[p][make_pair(trivial_left,trivial_left)] = prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.);
        for (const auto& iTerm: terms)
            this->add_term(iTerm);
    }


    /**
     * @brief Class constructor from a lattice and a vector of terms.
     * This method should be used for constructing operators != from the Hamiltonian.
     */
    TaggedMPOMaker(Lattice const& lat_, tag_vec const & i_, tag_vec const & f_,
                   std::shared_ptr<TagHandler<Matrix, SymmGroup> > th_, typename Model<Matrix, SymmGroup>::terms_type const& terms)
            : lat(lat_), identities(i_), fillings(f_), length(lat.size()), tag_handler(th_), prempo(length),
              trivial_left(prempo_key_type::trivial_left), trivial_right(prempo_key_type::trivial_right), leftmost_right(length),
              rightmost_left(0), finalized(false), verbose(false), core_energy(0.)
    {
        //for (size_t p = 0; p < length-1; ++p)
        //    prempo[p][make_pair(trivial_left,trivial_left)] = prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.);
        for (const auto& iTerm: terms)
            this->add_term(iTerm);
    }

    /** @brief Method to add a single term to the prempo object */
    void add_term(term_descriptor term)
    {
        std::sort(term.begin(), term.end(), pos_tag_lt());
        index_type nops = term.size();
        switch (nops) {
            case 1:
                add_1term(term);
                break;
            case 2:
                add_2term(term);
                break;
            default:
                add_generic_term(term);
                //add_nterm(term); /// here filling has to be done manually
                break;
        }
        leftmost_right = std::min(leftmost_right, term.rbegin()->first);
        rightmost_left = std::max(rightmost_left, term.begin()->first);
    }

    /** @bref Creates the MPO based on the tagged MPO object */
    MPO<Matrix, SymmGroup> create_mpo()
    {
        if (!finalized)
            finalize();
        MPO<Matrix, SymmGroup> mpo;
        mpo.reserve(length);
        typedef std::map<prempo_key_type, index_type> index_map;
        typedef typename index_map::iterator index_iterator;
        index_map left;
        left[trivial_left] = 0;
        typedef SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type> spin_desc_t;
        std::vector<spin_desc_t> left_spins(1);
        std::vector<index_type> LeftHerm(1);
        std::vector<int> LeftPhase(1,1);
        // Main loop over sites
        for (pos_t p = 0; p < length; ++p) {
            std::vector<tag_block> pre_tensor; pre_tensor.reserve(prempo[p].size());
            std::map<prempo_key_type, prempo_key_type> HermKeyPairs;
            std::map<prempo_key_type, std::pair<int,int> > HermitianPhases;
            index_map right;
            index_type r = 2;
            for (typename prempo_map_type::const_iterator it = prempo[p].begin(); it != prempo[p].end(); ++it) {
                // Remember that tag_block identifies a tuple with (size, size, tag_type and scaling factor)
                // The first two sizes should identify the operators acting on the left and right, respectively,
                // the tag_type is the operator and the scaling factor is the coefficient appearing in the definition
                // of the Hamiltonian
                prempo_key_type const& k1 = it->first.first;
                prempo_key_type const& k2 = it->first.second;
                prempo_value_type const& val = it->second;
                // Looks for the tag inside the left dictionary. Here an error is raised if the operator is not
                // found because, at each cycle, the left is filled with the right at the previous iteration
                index_iterator ll = left.find(k1);
                if (ll == left.end())
                    throw std::runtime_error("k1 not found!");
                // Looks for the tag in the right dictionary. If it has been not found,
                // updates the right dictionary
                index_iterator rr = right.find(k2);
                if (k2 == trivial_left && rr == right.end())
                    boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, 0) );
                else if (k2 == trivial_right && rr == right.end())
                    boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, 1) );
                else if (rr == right.end())
                    boost::tie(rr, boost::tuples::ignore) = right.insert( make_pair(k2, r++) );
                // Finalization
                index_type rr_dim = (p == length-1) ? 0 : rr->second;
                pre_tensor.push_back( tag_block(ll->second, rr_dim, val.first, val.second) );
                std::pair<int, int> phase;
                prempo_key_type ck2;
                boost::tie(ck2, phase) = conjugate_key(k2, p);
                if (!(k2 == ck2)){
                    HermKeyPairs[k2] = ck2;
                    HermitianPhases[k2] = phase;
                }
            }
            // Loads the dimensions of the left and right tags
            std::pair<index_type, index_type> rcd = rcdim(pre_tensor);
            // Spin-related part
            std::vector<spin_desc_t> right_spins(rcd.second);
            for (typename std::vector<tag_block>::const_iterator it = pre_tensor.begin(); it != pre_tensor.end(); ++it)
            {
                spin_desc_t out_spin = couple(left_spins[boost::tuples::get<0>(*it)],
                                              tag_handler->get_op(boost::tuples::get<2>(*it)).spin());
                index_type out_index = boost::tuples::get<1>(*it);
                assert(right_spins[out_index].get() == 0 || right_spins[out_index].get() == out_spin.get());
                right_spins[out_index] = out_spin;
            }
            // Locates the hermitian conjugate pairs
            std::vector<index_type> RightHerm(rcd.second);
            std::vector<int> RightPhase(rcd.second, 1);
            index_type cnt = 0;
            // For starting, everything is its own hermitian conjugate
            std::iota(RightHerm.begin(), RightHerm.end(), 0);
            for (typename std::map<prempo_key_type, prempo_key_type>::const_iterator
                            h_it = HermKeyPairs.begin(); h_it != HermKeyPairs.end(); ++h_it)
            {
                index_type romeo = right[h_it->first];
                index_type julia = right[h_it->second];
                if (romeo < julia)
                {
                    cnt++;
                    std::swap(RightHerm[romeo], RightHerm[julia]);
                    RightPhase[romeo] = HermitianPhases[h_it->first].first;
                    RightPhase[julia] = HermitianPhases[h_it->first].second;
                }
            }
            MPOTensor_detail::Hermitian h_(LeftHerm, RightHerm, LeftPhase, RightPhase);
            // Construction of the MPO tensor
            if (p == 0)
                mpo.push_back( MPOTensor<Matrix, SymmGroup>(1, rcd.second, pre_tensor,
                                 tag_handler->get_operator_table(), h_, left_spins, right_spins));
            else if (p == length - 1)
                mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, 1, pre_tensor,
                                 tag_handler->get_operator_table(), h_, left_spins, right_spins));
            else
                mpo.push_back( MPOTensor<Matrix, SymmGroup>(rcd.first, rcd.second, pre_tensor,
                                 tag_handler->get_operator_table(), h_, left_spins, right_spins));
            swap(left, right);
            swap(left_spins, right_spins);
            swap(LeftHerm, RightHerm);
            swap(LeftPhase, RightPhase);
            // Final Print
            if (verbose)
                maquis::cout << "MPO Bond " << p << ": " << rcd.second << "/" << cnt << std::endl;

        }
        mpo.setCoreEnergy(core_energy);
        return mpo;
    }

private:

    /**
     * @brief Adds a term composed by 1 SQ operators
     *
     * Note that if the SQ operator is the identity, the code just updates the core energy
     */
    void add_1term(term_descriptor const& term)
    {
        assert(term.size() == 1);
        /// Due to numerical instability: treat the core energy separately
        if (term.operator_tag(0) == identities[lat.get_prop<int>("type", term.position(0))]) {
            core_energy += std::real(term.coeff);
        }
        else {
            /// retrieve the actual operator from the tag table
            op_t current_op = tag_handler->get_op(term.operator_tag(0));
            current_op *= term.coeff;
            site_terms[term.position(0)] += current_op;
        }
    }

    /**
     * @brief Adds a term composed by 2 SQ operators
     * @param term_descriptor term: term to be added
     */
    void add_2term(term_descriptor const& term)
    {
        // Preliminary operations
        assert(term.size() == 2);
        auto Nt = lat.template get_prop<int>("NumTypes");
        SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
        // Number of fermionic operators of given type
        std::vector<int> v_nferm(Nt, 0);
        std::vector<bool> v_trivial_fill(Nt, true);  // This is used to keep track of the filling to be used per particle type
        std::vector<pos_t> v_pos(2);                 // Vector storing the position in which each operator acts
        std::vector<pos_t> v_part_type(2);           // Particle type for each site
        for (int i = 0; i < 2; ++i) {
          // Extract position and then type
          v_pos[i] = term.position(i);
          v_part_type[i] = lat.template get_prop<int>("ParticleType", std::vector<pos_t>{v_pos[i]});
          if (tag_handler->is_fermionic(term.operator_tag(i)))
            v_nferm[v_part_type[i]] += 1;
        }
        // First term
        prempo_key_type k1 = trivial_left;
        {
            int i = 0;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            prempo_key_type k2;
            k2.pos_op.push_back(term[i+1]);
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
              v_nferm[v_part_type[i]] -= 1;
            v_trivial_fill[v_part_type[i]] = (v_nferm[v_part_type[i]] % 2 == 0);
        }
        bool trivial_fill = !tag_handler->is_fermionic(term.operator_tag(1));
        insert_filling(term.position(0)+1, term.position(1), k1, v_trivial_fill, mpo_spin.get() > 1);
        // Second term
        {
            int i = 1;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            prempo_key_type k2 = trivial_right;
            insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
        }
        assert(mpo_spin.get() == 0); // H is a spin 0 operator
    }

    /**
     * @brief Adds a term composed by 3 SQ operators
     * @param term_descriptor term: term to be added
     */
    /*
    void add_3term(term_descriptor const& term)
    {
        assert(term.size() == 3);
        int nops = term.size();
        /// number of fermionic operators
        int nferm = 0;
        for (int i = 0; i < nops; ++i) {
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm += 1;
        }
        SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
        prempo_key_type k1 = trivial_left;
        std::vector<pos_op_type> ops_left;
        /// op_0
        {
            int i = 0;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            prempo_key_type k2;
            k2.pos_op.push_back(to_pair(term[i])); // k2: applied operator
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm -= 1;
            bool trivial_fill = (nferm % 2 == 0);
            insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
        }
        /// op_1
        {
            int i = 1;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            prempo_key_type k2;
            k2.pos_op.push_back(to_pair(term[i+1])); // k2: future operators
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), maquis::real(term.coeff)), detach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm -= 1;
            bool trivial_fill = (nferm % 2 == 0);
            insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
        }
        /// op_2
        {
            int i = 2;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            insert_operator(term.position(i), make_pair(k1, trivial_right), prempo_value_type(term.operator_tag(i), 1.), attach);
        }
        assert(mpo_spin.get() == 0); // H is a spin 0 operator
    }
    */

    /*
    void add_4term(term_descriptor const& term)
    {
        assert(term.size() == 4);
        int nops = term.size();
        /// number of fermionic operators
        int nferm = 0;
        for (int i = 0; i < nops; ++i) {
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm += 1;
        }
        SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
        prempo_key_type k1 = trivial_left;
        std::vector<pos_op_type> ops_left;
        /// op_0, op_1
        for (int i = 0; i < 2; ++i) {
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            ops_left.push_back(to_pair(term[i])); prempo_key_type k2(ops_left);
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm -= 1;
            bool trivial_fill = (nferm % 2 == 0);
            insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
        }
        /// op_2
        {
            int i = 2;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            prempo_key_type k2;
            k2.pos_op.push_back(to_pair(term[3]));
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), maquis::real(term.coeff)), detach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                nferm -= 1;
            bool trivial_fill = (nferm % 2 == 0);
            insert_filling(term.position(i)+1, term.position(i+1), k1, trivial_fill, (mpo_spin.get() > 1) ? term.full_identity : -1);
        }
        /// op_3
        {
            int i = 3;
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            insert_operator(term.position(i), make_pair(k1, trivial_right), prempo_value_type(term.operator_tag(i), 1.), attach);
        }
        assert(mpo_spin.get() == 0); // H is a spin 0 operator
    }
    */

    /**
     * @brief Add a general n-body term using the optimized MPO construction.
     *
     * This method generalizes the [add_3body] and [add_4body] operator to arbitrarily
     * complex many-body coupling terms
     * When the boolean prefer_fork is true, if the number of operators is odd uses
     * one additional fork label, otherwise uses the merge one.
     */
    void add_generic_term(const term_descriptor& term, bool const& prefer_fork=true)
    {
        std::size_t nops = term.size();
        // Case nops == 2 is managed separately, inside [add_2term]
        assert ( nops != 2 ) ;
        auto Nt = lat.template get_prop<int>("NumTypes");
        // Number of fermionic operators of given type
        std::vector<int> v_nferm(Nt, 0);
        std::vector<bool> v_trivial_fill(Nt, true);
        std::vector<pos_t> v_pos(nops);
        std::vector<pos_t> v_part_type(nops);
        for (int i = 0; i < nops; ++i) {
            // Extract position and then type
            v_pos[i] = term.position(i);
            v_part_type[i] = lat.template get_prop<int>("ParticleType", std::vector<pos_t>{v_pos[i]});
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                v_nferm[v_part_type[i]] += 1;
        }
        SpinDescriptor<typename symm_traits::SymmType<SymmGroup>::type > mpo_spin;
        std::vector<pos_op_type> ops_left, ops_right;
        // Thresh is the index after which we switch from fork to merge behaviour
        int thresh;
        if (nops % 2 == 1 || prefer_fork)
            thresh = nops/2;
        else
            thresh = nops/2 - 1;
        // == FORKING OPERATORS ==
        prempo_key_type k1 = trivial_left;
        for (std::size_t i = 0; i < thresh; ++i) {
          mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
          ops_left.push_back(term[i]);
          prempo_key_type k2(ops_left);
          k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), attach);
          // Checks how many fermionic operators are left - if the number is odd, will insert the filling,
          // otherwise will insert the identity
          if (tag_handler->is_fermionic(term.operator_tag(i)))
            v_nferm[v_part_type[i]] -= 1;
          v_trivial_fill[v_part_type[i]] = (v_nferm[v_part_type[i]] % 2 == 0);
          // -> if types at  term.position(i) and  term.position(i+1) are different insert trivial fill
          insert_filling(term.position(i)+1, term.position(i+1), k1, v_trivial_fill, mpo_spin.get() > 1);
        }
        // == MIDDLE OPERATOR ==
        // Note that in this case we use the detach modality
        prempo_key_type k2;
        for (std::size_t j = thresh+1; j < nops; j++)
            ops_right.push_back(term[j]);
        k2 = prempo_key_type(ops_right);
        mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(thresh))).spin());
        k1 = insert_operator(term.position(thresh), make_pair(k1, k2), prempo_value_type(term.operator_tag(thresh), term.coeff), detach);
        // Extract position and then type
        if (tag_handler->is_fermionic(term.operator_tag(thresh)))
            v_nferm[v_part_type[thresh]] -= 1;
        v_trivial_fill[v_part_type[thresh]] = (v_nferm[v_part_type[thresh]] % 2 == 0);
        insert_filling(term.position(thresh)+1, term.position(thresh+1), k1, v_trivial_fill, mpo_spin.get() > 1);
        // == MERGE OPERATOR ==
        for (std::size_t i = thresh+1; i < nops; i++) {
            // Extract position and then type
            ops_right.resize(0);
            prempo_key_type k2;
            if ( i == nops-1 ) {
                k2 = trivial_right;
            } else {
                for (int j = i+1; j < nops; j++)
                    ops_right.push_back(term[j]);
                k2 = prempo_key_type(ops_right);
            }
            mpo_spin = couple(mpo_spin, (tag_handler->get_op(term.operator_tag(i))).spin());
            k1 = insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.),
                                 attach);
            if (tag_handler->is_fermionic(term.operator_tag(i)))
                v_nferm[v_part_type[i]] -= 1;
            if ( i != nops-1 ) {
                v_trivial_fill[v_part_type[i]] = (v_nferm[v_part_type[i]] % 2 == 0);
                insert_filling(term.position(i)+1, term.position(i+1), k1, v_trivial_fill, mpo_spin.get() > 1);
            }
        }
        // Final term, the only one where the coefficient is actually employed
        assert(mpo_spin.get() == 0); // H is a spin 0 operator
    }

    /**
     * @brief Adds a term to the MPO using the naive construction
     *
     * Note that in this case the filling operators must be included in the terms_,
     * since they are not added by the `tagged_mpo_maker`.
     */
    void add_nterm(term_descriptor const& term)
    {
        int nops = term.size();
        assert( nops > 2 );
        static index_type next_offset = 0;
        index_type current_offset = (next_offset++);
        prempo_key_type k1 = trivial_left;
        prempo_key_type k2(prempo_key_type::bulk_no_merge, current_offset);
        k2.pos_op.push_back(term[nops-1]);
        {
            int i = 0;
            insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), term.coeff), detach);
            k1 = k2;
            if (i < nops-1 && term.position(i)+1 != term.position(i+1))
                throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
        }
        //
        for (int i = 1; i < nops; ++i) {
            if (i == nops-1)
                k2 = trivial_right;
            //
            insert_operator(term.position(i), make_pair(k1, k2), prempo_value_type(term.operator_tag(i), 1.), detach);
            if (i < nops-1 && term.position(i)+1 != term.position(i+1))
                throw std::runtime_error("for n > 4 operators filling is assumed to be done manually. the list of operators contains empty sites.");
        }
    }

    /**
     * @brief Insert proper the filling operators between sites i and j.
     * @param i starting site
     * @param j ending site
     * @param k prempo_key_type associated with the operator to be propagated
     * @param trivial_fill vector of bool with size == number of particle types indicating whether the identity or
     * the filling operator should be used
     * @param isSpinLargerThanOne if true, uses the identity full operator instead of the "simpler" identity
     */
    void insert_filling(pos_t i, pos_t j, prempo_key_type k, std::vector<bool> trivial_fill, bool isSpinLargerThanOne)
    {
        for (; i < j; ++i) {
            auto typei = lat.get_prop<int>("type", i);
            auto particleTypei = lat.get_prop<int>("ParticleType", i);
            tag_type use_ident = (isSpinLargerThanOne) ? identities_full[typei] : identities[typei];
            tag_type op = (trivial_fill[particleTypei]) ? use_ident : fillings[typei];
            //std::pair<typename prempo_map_type::iterator,bool> ret = prempo[i].insert( make_pair(make_pair(k,k), prempo_value_type(op, 1.)) );
            //if (!ret.second && ret.first->second.first != op)
            if (prempo[i].count(make_pair(k,k)) == 0) {
                auto ret = prempo[i].insert( make_pair(make_pair(k,k), prempo_value_type(op, 1.)) );
            }
            else {
                if (prempo[i].find(make_pair(k,k))->second != prempo_value_type(op, 1.))
                throw std::runtime_error("Pre-existing term at site "+std::to_string(i)+ ". Needed "+std::to_string(op)
                                            + ", found "+std::to_string(prempo[i].find(make_pair(k,k))->second.first));
            }
        }
    }

    /**
     * @brief Method to insert an element in the prempo objecj
     * @param int p: site on which the operator acts.
     * @param pair<prempo_key_type, prempo_key_type> kk: element to add
     * @param prempo_value_type val: element to be added on site p
     * @param merge_kind merge_behaviour: detach if a new branch should not be created, attach otherwise
     */
    prempo_key_type insert_operator(pos_t p, std::pair<prempo_key_type, prempo_key_type> kk, prempo_value_type val,
                                    merge_kind merge_behavior=detach)
    {
        /// merge_behavior == detach: a new branch will be created, in case op already exist, an offset is used
        /// merge_behavior == attach: if operator tags match, keep the same branch
        if (merge_behavior == detach)
            prempo[p].insert( make_pair(kk, val) );
        else
            if (prempo[p].count(kk) == 0)
                prempo[p].insert( make_pair(kk, val) );
        return kk.second;
    }

    /**
     * @brief Performs final operations on the prempo object.
     *
     * Manages the one-site operators and add identity at the two extreme of the lattice (if it was not
     * done before)
     */
    void finalize()
    {
        /// site terms
        std::pair<prempo_key_type,prempo_key_type> kk = make_pair(trivial_left,trivial_right);
        for (typename std::map<pos_t, op_t>::const_iterator it = site_terms.begin();
             it != site_terms.end(); ++it) {
            tag_type site_tag = tag_handler->register_op(it->second, tag_detail::bosonic);
            //std::pair<typename prempo_map_type::iterator,bool> ret;
            //ret = prempo[it->first].insert( make_pair( kk, prempo_value_type(site_tag,1.) ) );
            typename prempo_map_type::iterator ret;
            ret = prempo[it->first].insert( make_pair( kk, prempo_value_type(site_tag, 1.) ) );
            if (prempo[it->first].count(ret->first) != 1)
                throw std::runtime_error("another site term already existing!");
        }
        // fill with ident from the begin
        for (size_t p = 0; p < rightmost_left; ++p)
            prempo[p].insert(make_pair(make_pair(trivial_left,trivial_left),
                                       prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.)));
        /// fill with ident until the end
        for (size_t p = leftmost_right+1; p < length; ++p)
            prempo[p].insert(make_pair(make_pair(trivial_right,trivial_right),
                                       prempo_value_type(identities[lat.get_prop<int>("type",p)], 1.)));
        finalized = true;
    }

    /**
     * @brief Finds the conjugate key of a given entry.
     *
     * Note that, in the model construction, we register pairs of Hermitian conjugate *elementary*
     * operators (i.e., operators that act on a single site).
     */
    std::pair<prempo_key_type, std::pair<int, int> > conjugate_key(prempo_key_type k, pos_t p)
    {
        // Defines a function to get the number of particles of a given symmetry charge
        typename SymmGroup::subcharge (*np)(typename SymmGroup::charge) = &SymmGroup::particleNumber;
        //if (k.pos_op.size() > 1)
        //    return std::make_pair(k, std::make_pair(1,1));
        prempo_key_type conj = k;
        // Loop over all the elements of the prempo_key_type operator and finds the corresponding
        // Hermitian conjugate. Note that it is sufficient that one of the key has no corresponding
        // conjugate to not register the hermitian conjugate.
        for (tag_type i = 0; i < k.pos_op.size(); ++i) {
            //if (k.pos_op[i].second == tag_handler->herm_conj(k.pos_op[i].second))
            //    return std::make_pair(k, std::make_pair(1,1));
            conj.pos_op[i].second = tag_handler->herm_conj(k.pos_op[i].second);
        }
        // Calculates the phase correction for the hermitian conjugate of SU(2) operators.
        std::pair<int, int> phase(1,1);
        if ( k.pos_op.size() == 1) {
            // merge type operator ahead of current position
            if ( p < k.pos_op[0].first ) {
                SiteOperator<Matrix, SymmGroup> const & op1 = tag_handler->get_op(k.pos_op[0].second);
                typename SymmGroup::subcharge pdiff = np(op1.basis().left_charge(0)) - np(op1.basis().right_charge(0));
                if ( pdiff == 1) //  creator
                    phase = std::make_pair(1, -1);
                else if ( pdiff == -1) // destructor
                    phase = std::make_pair(-1, 1);
            }
            else {
                SiteOperator<Matrix, SymmGroup> const & op1 = tag_handler->get_op(k.pos_op[0].second);
                if ( op1.spin().get() == 1) // creator or destructor
                    phase = std::make_pair(-1, 1);
            }
        }
        //
        if ( k.pos_op.size() == 2) {
            SiteOperator<Matrix, SymmGroup> const & op1 = tag_handler->get_op(k.pos_op[0].second);
            SiteOperator<Matrix, SymmGroup> const & op2 = tag_handler->get_op(k.pos_op[1].second);
            // if k contains (c^dag c)_S=0 or (c c^dag)_S=0
            if (op1.spin().get() == 1 && op2.spin().get() == 1 && op2.spin().action() == -1
                && np(op1.basis().left_charge(0)) - np(op1.basis().right_charge(0)) ==
                - (np(op2.basis().left_charge(0)) - np(op2.basis().right_charge(0)))
               )
                phase = std::make_pair(-1,-1);
            // if k contains (c^dag c^dag)_S=1 or (c c)_S=1
            if (op1.spin().get() == 1 && op2.spin().get() == 1 && op2.spin().action() == 1
                && np(op1.basis().left_charge(0)) - np(op1.basis().right_charge(0)) ==
                  (np(op2.basis().left_charge(0)) - np(op2.basis().right_charge(0)))
               )
                phase = std::make_pair(-1,-1);
        }
        return std::make_pair(conj, phase);
    }

private:
    Lattice const& lat;
    tag_vec identities, identities_full, fillings;
    pos_t length;
    std::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    std::vector<prempo_map_type> prempo;
    prempo_key_type trivial_left, trivial_right;
    std::map<pos_t, op_t> site_terms;
    pos_t leftmost_right, rightmost_left;
    bool finalized, verbose;
    double core_energy;
};

} // namespace generate_mpo

#endif
