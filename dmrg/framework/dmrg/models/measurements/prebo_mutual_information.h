/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_PREBO_MUTUAL_INFORMATION_H
#define MAQUIS_DMRG_PREBO_MUTUAL_INFORMATION_H

/* internal */
#include "dmrg/models/measurement.h"
#include "dmrg/utils/checks.h"
#include "measurements_details.h"
#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/prebo/prebo_TermGenerator.hpp"
/* external */
#include <cmath>
#include <alps/numeric/matrix/algorithms.hpp>

#ifdef DMRG_PREBO

namespace measurements {

    /**
     * @brief Measurement class associated with the PreBO Mutual information
     * @tparam Matrix numeric matrix class
     * @tparam N integer associated with the NU1 symmetry class
     */
    template <class Matrix, int N>
    class PreBOMutualInformation : public measurement<Matrix, NU1_template<N>> {
    public:
        // Types definition
        using NU1 = NU1_template<N>;
        using tag_type = typename generate_mpo::OperatorTagTerm<Matrix, NU1_template<N>>::tag_type;
        using tag_vec = std::vector<tag_type>;
        using base = measurement<Matrix, NU1>;

    public:
        /**
         * @brief Class constructor
         * @param parms Parameter container
         * @param lat_ DMRG lattice
         * @param identities vector with the tags associated with the identity operators
         * @param fillings vector with the tags associated with the filling operators
         * @param ptr_term_generator_ tag generator object (required by the MPO management in the NU1 class)
         */
        PreBOMutualInformation(BaseParameters& parms, Lattice const& lat_, const tag_vec& identities, const tag_vec& fillings,
                         std::shared_ptr<prebo::TermGenerator<Matrix, N>> ptr_term_generator_, bool verbose=true)
            : base("mutinf"), parms(parms), lat(lat_), identities(identities), fillings(fillings),
              ptr_term_generator(ptr_term_generator_), verbose_(verbose)
        { }

        /** @brief Class destructor */
        ~PreBOMutualInformation()=default;

        /**
         * @brief Implementation of virtual method in base class.
         * @param ket_mps input MPS.
         * @param rmps reduced MPS (required by the virtual class interface).
         */
        void evaluate(MPS<Matrix, NU1> const& ket_mps, boost::optional<reduced_mps<Matrix, NU1> const&> rmps = boost::none) {

            this->vector_results.clear();
            this->labels.clear();
            this->labels_num.clear();
            std::vector<typename MPS<Matrix, NU1>::scalar_type> dct;
            std::vector<std::vector<Lattice::pos_t>> num_labels;

            auto vec_orbitals = lat.template get_prop<std::vector<int>>("vec_orbitals");
            auto isFermion = lat.template get_prop<std::vector<bool>>("isFermion");
            auto inv_order = lat.template get_prop<std::vector<Lattice::pos_t>>("inv_order");
            auto L = std::accumulate(vec_orbitals.begin(), vec_orbitals.end(), 0);
            Matrix entropyTable(L, L, 0.0);
            // Loop over all sites.
            if (verbose_) {
                std::cout << "\nCalulating orbital entropies:\n\n";
                std::cout << "         i-mu         j-nu      Tr            S\n";
                std::cout << "----------------------------------------------------------\n";
            }
            for (size_t i=0; i<vec_orbitals.size(); i++) {
                for (size_t j = 0; j < vec_orbitals.size(); j++) {
                    for (size_t mu = 0; mu < vec_orbitals.at(i); mu++) {
                        for (size_t nu = 0; nu < vec_orbitals.at(j); nu++) {

                            if (i == j && nu > mu) continue;
                            if (j > i) continue;

                            double trace = 0;
                            auto site1 = NBodyTerm::retrieve_abs_index(mu, i, vec_orbitals, inv_order);
                            auto site2 = NBodyTerm::retrieve_abs_index(nu, j, vec_orbitals, inv_order);
                            if (site2 > site1) {
                                site1 = NBodyTerm::retrieve_abs_index(nu, j, vec_orbitals, inv_order);
                                site2 = NBodyTerm::retrieve_abs_index(mu, i, vec_orbitals, inv_order);
                            }

                            // **********
                            //   1o-RDM
                            // **********
                            if (i == j && isFermion[i] && mu == nu) {
                                std::vector<double> OneOrbRDM(4, 0.);
                                // rho = diag(O_1, O_6, O_11, O_16)
                                std::vector<size_t> vals{1, 6, 11, 16};
                                for (size_t iter = 0; iter < 4; iter++) {
                                    auto terms = ptr_term_generator->generate_termsTransitionOp(i, site1,
                                                                                                vals.at(iter));
                                    generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                   ptr_term_generator->getTagHandler(),
                                                                                   terms);
                                    MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                    MPO<Matrix, NU1> mpoc = mpo;
                                    OneOrbRDM.at(iter) = maquis::real(expval(ket_mps, mpoc));
                                }
                                double entropy = 0;
                                for (auto const w: OneOrbRDM) {
                                    if (!std::isnan(w))
                                        trace += w;
                                    if (w > 1e-15)
                                        entropy += w * std::log(w);
                                }
                                // Finalize:
                                entropyTable(site1, site1) = -entropy;
                            }
                                // **********
                                //   2o-RDM
                                // **********
                                // First case: same particle type and fermionic
                                // The orbital entropy is of course symmetric with permuting mu and nu.
                            else if (i == j && isFermion[i]) {//} && mu>nu) {
                                // In this case we have 9 blocks which are
                                //      N      S    dim
                                // 1.   0      0      1
                                // 2.   1   -1/2      2
                                // 3.   1    1/2      2
                                // 4.   2     -1      1
                                // 5.   2      0      4
                                // 6.   2      1      1
                                // 7.   3   -1/2      2
                                // 8.   3    1/2      2
                                // 4.   4      0      1
                                double entropy = 0;
                                // Block 1
                                // N=0 S=0
                                // 1/1
                                {
                                    double w = 0;
                                    auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2, 1,
                                                                                                1);
                                    generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                   ptr_term_generator->getTagHandler(),
                                                                                   terms);
                                    MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                    MPO<Matrix, NU1> mpoc = mpo;
                                    w = maquis::real(expval(ket_mps, mpoc));
                                    if (!std::isnan(w))
                                        trace += w;
                                    if (w > 1e-15)
                                        entropy += w * std::log(w);
                                }
                                // Block 2
                                // N=1 S=-1/2
                                // 1/6 2/5
                                //     6/1
                                {
                                    Matrix block(2, 2, 0.0);
                                    // 1/6
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    1, 6);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 0) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 6/1
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    6, 1);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 1) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 2/5 and 5/2
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    2, 5);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 1) = maquis::real(expval(ket_mps, mpoc));
                                        block(1, 0) = block(0, 1);
                                    }
                                    Matrix evecs(block.num_rows(), block.num_rows(), 0.0);
                                    alps::numeric::vector<double> evals(block.num_rows(), 0.0);
                                    alps::numeric::syev(block, evecs, evals);
                                    for (size_t it = 0; it < evals.size(); it++) {
                                        if (!std::isnan(evals(it)))
                                            trace += evals(it);
                                        if (evals(it) > 1e-15)
                                            entropy += evals(it) * std::log(evals(it));
                                    }
                                } // block2
                                // Block 3
                                // N=1 S=1/2
                                // 1/11   3/9
                                //       11/1
                                {
                                    Matrix block(2, 2, 0.0);
                                    // 1/11
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    1, 11);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 0) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 11/1
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    11, 1);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 1) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 3/9 and 9/3
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    3, 9);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 1) = maquis::real(expval(ket_mps, mpoc));
                                        block(1, 0) = block(0, 1);
                                    }
                                    Matrix evecs(block.num_rows(), block.num_rows(), 0.0);
                                    alps::numeric::vector<double> evals(block.num_rows(), 0.0);
                                    alps::numeric::syev(block, evecs, evals);
                                    for (size_t it = 0; it < evals.size(); it++) {
                                        if (!std::isnan(evals(it)))
                                            trace += evals(it);
                                        if (evals(it) > 1e-15)
                                            entropy += evals(it) * std::log(evals(it));
                                    }
                                } //block3
                                // Block 4
                                // N=2 S=-1
                                // 6/6
                                {
                                    double w = 0;
                                    auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2, 6,
                                                                                                6);
                                    generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                   ptr_term_generator->getTagHandler(),
                                                                                   terms);
                                    MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                    MPO<Matrix, NU1> mpoc = mpo;
                                    w = maquis::real(expval(ket_mps, mpoc));
                                    if (!std::isnan(w))
                                        trace += w;
                                    if (w > 1e-15)
                                        entropy += w * std::log(w);
                                }
                                // Block 5
                                // N=2 S=0
                                // 1/16   2/15    3/14    4/13
                                //        6/11    7/10    8/9
                                //                11/6    12/5
                                //                        16/1
                                {
                                    Matrix block(4, 4, 0.0);
                                    // Diagonal elements:
                                    // 0,0: 1/16
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    1, 16);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 0) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 1,1: 6/11
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    6, 11);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 1) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 2,2: 11/6
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    11, 6);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(2, 2) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 3,3: 16/1
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    16, 1);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(3, 3) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // Off-diagonal
                                    // 0,1 and 1,0: 2/15
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    2, 15);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 1) = maquis::real(expval(ket_mps, mpoc));
                                        block(1, 0) = block(0, 1);
                                    }
                                    // 0,2 and 2,0: 3/14
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    3, 14);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 2) = maquis::real(expval(ket_mps, mpoc));
                                        block(2, 0) = block(0, 2);
                                    }
                                    // 0,3 and 3,0: 4/13
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    4, 13);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 3) = maquis::real(expval(ket_mps, mpoc));
                                        block(3, 0) = block(0, 3);
                                    }
                                    // 1,2 and 2,1: 7/10
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    7, 10);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 2) = maquis::real(expval(ket_mps, mpoc));
                                        block(2, 1) = block(1, 2);
                                    }
                                    // 1,3 and 3,1: 8/9
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    8, 9);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 3) = maquis::real(expval(ket_mps, mpoc));
                                        block(3, 1) = block(1, 3);
                                    }
                                    // 2,3 and 3,2: 12/5
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    12, 5);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(2, 3) = maquis::real(expval(ket_mps, mpoc));
                                        block(3, 2) = block(2, 3);
                                    }
                                    Matrix evecs(block.num_rows(), block.num_rows(), 0.0);
                                    alps::numeric::vector<double> evals(block.num_rows(), 0.0);
                                    alps::numeric::syev(block, evecs, evals);
                                    for (size_t it = 0; it < evals.size(); it++) {
                                        if (!std::isnan(evals(it)))
                                            trace += evals(it);
                                        if (evals(it) > 1e-15)
                                            entropy += evals(it) * std::log(evals(it));
                                    }
                                } //block5
                                // Block 6
                                // N=2 S=1
                                // 11/11
                                {
                                    double w = 0;
                                    auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2, 11,
                                                                                                11);
                                    generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                   ptr_term_generator->getTagHandler(),
                                                                                   terms);
                                    MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                    MPO<Matrix, NU1> mpoc = mpo;
                                    w = maquis::real(expval(ket_mps, mpoc));
                                    if (!std::isnan(w))
                                        trace += w;
                                    if (w > 1e-15)
                                        entropy += w * std::log(w);
                                }
                                // Block 7
                                // N=3 S=-1/2
                                // 6/16   8/14
                                //        16/6
                                {
                                    Matrix block(2, 2, 0.0);
                                    // 6/16
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    6, 16);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 0) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 16/6
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    16, 6);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 1) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 8/14 and 14/8
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    8, 14);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 1) = maquis::real(expval(ket_mps, mpoc));
                                        block(1, 0) = block(0, 1);
                                    }
                                    Matrix evecs(block.num_rows(), block.num_rows(), 0.0);
                                    alps::numeric::vector<double> evals(block.num_rows(), 0.0);
                                    alps::numeric::syev(block, evecs, evals);
                                    for (size_t it = 0; it < evals.size(); it++) {
                                        if (!std::isnan(evals(it)))
                                            trace += evals(it);
                                        if (evals(it) > 1e-15)
                                            entropy += evals(it) * std::log(evals(it));
                                    }
                                } //block7
                                // Block 8
                                // N=3 S=1/2
                                // 11/16  12/15
                                //        16/11
                                {
                                    Matrix block(2, 2, 0.0);
                                    // 11/16
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    11, 16);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 0) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 16/11
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    16, 11);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(1, 1) = maquis::real(expval(ket_mps, mpoc));
                                    }
                                    // 12/15 and 15/12
                                    {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    12, 15);
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        block(0, 1) = maquis::real(expval(ket_mps, mpoc));
                                        block(1, 0) = block(0, 1);
                                    }
                                    Matrix evecs(block.num_rows(), block.num_rows(), 0.0);
                                    alps::numeric::vector<double> evals(block.num_rows(), 0.0);
                                    alps::numeric::syev(block, evecs, evals);
                                    for (size_t it = 0; it < evals.size(); it++) {
                                        if (!std::isnan(evals(it)))
                                            trace += evals(it);
                                        if (evals(it) > 1e-15)
                                            entropy += evals(it) * std::log(evals(it));
                                    }
                                } //block8
                                // Block 9
                                // N=4 S=0
                                // 16/16
                                {
                                    double w = 0;
                                    auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2, 16,
                                                                                                16);
                                    generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                   ptr_term_generator->getTagHandler(),
                                                                                   terms);
                                    MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                    MPO<Matrix, NU1> mpoc = mpo;
                                    w = maquis::real(expval(ket_mps, mpoc));
                                    if (!std::isnan(w))
                                        trace += w;
                                    if (w > 1e-15)
                                        entropy += w * std::log(w);
                                }
                                // Finalize:
                                entropyTable(site1, site2) = -entropy;
                            } else if (i != j && isFermion[i] && isFermion[j]) {
                                double entropy = 0;
                                std::vector<size_t> vals{1, 6, 11, 16};
                                for (size_t iter1 = 0; iter1 < 4; iter1++)
                                    for (size_t iter2 = 0; iter2 < 4; iter2++) {
                                        auto terms = ptr_term_generator->generate_termsTransitionOp(i, j, site1, site2,
                                                                                                    vals.at(iter1),
                                                                                                    vals.at(iter2));
                                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings,
                                                                                       ptr_term_generator->getTagHandler(),
                                                                                       terms);
                                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                                        MPO<Matrix, NU1> mpoc = mpo;
                                        auto w = maquis::real(expval(ket_mps, mpoc));
                                        if (!std::isnan(w))
                                            trace += w;
                                        if (w > 1e-15)
                                            entropy += w * std::log(w);
                                    }
                                // Finalize:
                                entropyTable(site1, site2) = -entropy;
                            }
                            if (verbose_) {
                                std::cout << std::setw(10) << i << "-" << mu
                                          << std::setw(10) << j << "-" << nu << std::setw(10) << trace
                                          << std::setw(20) << entropyTable(site1, site2) << "\n";
                            }
                        }
                    }
                }
            }
            // save results and labels
            for (int row=0; row < L; row++) {
                for (int col = 0; col <= row; col++) {
                    // Single-orbital entropies on diagonal amd mutual information on off-diagonal
                    if (row==col) {
                        dct.push_back(entropyTable(row,row));
                    }
                    else {
                        dct.push_back(0.5 * (entropyTable(row, row) + entropyTable(col, col) - entropyTable(row, col)));
                    }
                    num_labels.push_back({row, col});
                }
            }
            std::vector<std::string> lbt = label_strings(num_labels);
            {
                this->vector_results.reserve(this->vector_results.size() + dct.size());
                std::copy(dct.begin(), dct.end(), std::back_inserter(this->vector_results));

                this->labels.reserve(this->labels.size() + dct.size());
                std::copy(lbt.begin(), lbt.end(), std::back_inserter(this->labels));

                this->labels_num.reserve(this->labels_num.size() + dct.size());
                std::copy(num_labels.begin(), num_labels.end(), std::back_inserter(this->labels_num));
            }
            if (parms.is_set("mutinffile")) {
                std::string fname= parms["mutinffile"].str() + ".csv";
                std::cout << std::endl;
                std::cout << "Writing mutual information on disk ... " << std::endl;
                std::ofstream file(fname);
                for (auto row=0; row<L; ++row) {
                    for (auto col=0; col<L-1; ++col) {
                        file << entropyTable(row,col) << ",";
                    }
                    file << entropyTable(row,L-1) << ",\n";
                }
                file.close();
            }
        }

    protected:
        measurement<Matrix, NU1>* do_clone() const
        {
            return new PreBOMutualInformation(*this);
        }

    private:
        BaseParameters& parms;
        const Lattice& lat;
        tag_vec identities, fillings;
        std::shared_ptr<prebo::TermGenerator<Matrix, N>> ptr_term_generator;
        bool verbose_;
    };
}

#endif // DMRG_PREBO

#endif // MAQUIS_DMRG_PREBO_MUTUAL_INFORMATION_H
