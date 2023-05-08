/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_PREBO_PARTICLE_RDM_H
#define MAQUIS_DMRG_PREBO_PARTICLE_RDM_H

/* internal */
#include "dmrg/models/measurement.h"
#include "dmrg/utils/checks.h"
#include "measurements_details.h"
#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/models/prebo/nu1/model.hpp"
#include "dmrg/models/prebo/prebo_TermGenerator.hpp"
/* external */
#include <alps/numeric/matrix/algorithms.hpp>

#ifdef DMRG_PREBO

namespace measurements {

    /**
     * @class PreBOParticleRDM
     * @brief Implementation of the virtual base class for the calculation of the particle RDM with the pre-BO model.
     * @tparam Matrix Matrix type underlying the MPS
     * @tparam N integer associated with the NU1 Symm Group
     */
    template <class Matrix, int N>
    class PreBOParticleRDM : public measurement<Matrix, NU1_template<N>> {
    public:
        // Types definition
        using NU1 = NU1_template<N>;
        using tag_type = typename generate_mpo::OperatorTagTerm<Matrix, NU1>::tag_type;
        using tag_vec = std::vector<tag_type>;
        using base = measurement<Matrix, NU1>;

        /**
         * @brief Class constructor
         * @param parms parameter container
         * @param lat_ lattice
         * @param identities tags associated with the identity operators
         * @param fillings tags associated with the filling operators
         * @param ptr_term_generator_ term generator object (see prebo_TermGenerator.hpp)
         */
        PreBOParticleRDM(BaseParameters& parms, Lattice const& lat_, const tag_vec& identities, const tag_vec& fillings,
                         std::shared_ptr<prebo::TermGenerator<Matrix, N>> ptr_term_generator_) 
            : base("oneptdm"), parms(parms), lat(lat_), identities(identities), fillings(fillings), ptr_term_generator(ptr_term_generator_) 
        { }

        /** @brief Default destructor */
        ~PreBOParticleRDM()=default;

        /**
         * @brief Implementation of virtual method [evaluate] required by base class.
         * @param ket_mps Input MPS
         * @param rmps reduced MPS (needed to match the interface)
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
            for (int nt=0; nt<vec_orbitals.size(); ++nt) {
                size_t dim = vec_orbitals.at(nt);
                Matrix rdm(vec_orbitals.at(nt), vec_orbitals.at(nt), 0.0);
                for (int mu = 0; mu < dim; mu++) {
                    for (int nu = 0; nu <= mu; nu++) {
                        auto terms = ptr_term_generator->generate_terms1RDM(nt, mu, nu);
                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings, ptr_term_generator->getTagHandler(), terms);
                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                        MPO<Matrix, NU1> mpoc = mpo;
                        rdm(mu, nu) = maquis::real(expval(ket_mps, mpoc));
                        if (nu != mu)
                            rdm(nu, mu) = rdm(mu, nu);
                        num_labels.push_back({nt, mu, nu});
                        dct.push_back(rdm(mu,nu));
                    }
                }
                std::vector<std::string> lbt = label_strings(num_labels);
                // save results and labels
                {
                    this->vector_results.reserve(this->vector_results.size() + dct.size());
                    std::copy(dct.begin(), dct.end(), std::back_inserter(this->vector_results));

                    this->labels.reserve(this->labels.size() + dct.size());
                    std::copy(lbt.begin(), lbt.end(), std::back_inserter(this->labels));

                    this->labels_num.reserve(this->labels_num.size() + dct.size());
                    std::copy(num_labels.begin(), num_labels.end(), std::back_inserter(this->labels_num));
                }

                #ifndef NDEBUG
                    std::cout << std::endl;
                    std::cout << "The one-body RDM is" << std::endl;
                    std::cout << rdm << std::endl;
                    std::cout << std::endl;
                #endif
                Matrix evecs(vec_orbitals.at(nt), vec_orbitals.at(nt), 0.0);
                alps::numeric::vector<double> evals(vec_orbitals.at(nt), 0.0);
                alps::numeric::syev(rdm, evecs, evals);
                std::cout << "Eigenvalues of one-body RDM of type " << nt << " are:" << std::endl;
                for (auto i=0; i< vec_orbitals.at(nt); ++i) {
                    std::cout << std::setw(16) << evals(i) << std::endl;
                }
                std::cout << "The trace is:      " << alps::numeric::trace(rdm) << std::endl;
                std::cout << "The norm is:    " << norm(ket_mps) << std::endl;
                if (parms.is_set("rdmfile")) {
                    std::string fname= parms["rdmfile"].str() + "_" + std::to_string(nt) + ".csv";
                    std::cout << "Writing 1-body rdm on disk ... " << std::endl;
                    std::ofstream file(fname);
                    for (auto row=0; row<vec_orbitals.at(nt); ++row) {
                        for (auto col=0; col<vec_orbitals.at(nt)-1; ++col) {
                            file << rdm(row,col) << ",";
                        }
                        file << rdm(row,vec_orbitals.at(nt)-1) << "\n";
                    }
                    file.close();
                }
            }
        }

    protected:
        measurement<Matrix, NU1>* do_clone() const
        {
            return new PreBOParticleRDM(*this);
        }

    private:
        BaseParameters& parms;
        const Lattice& lat;
        tag_vec identities, fillings;
        std::shared_ptr<prebo::TermGenerator<Matrix, N>> ptr_term_generator;
    };
}

#endif // DMRG_PREBO

#endif // MAQUIS_DMRG_PREBO_PARTICLE_RDM_H
