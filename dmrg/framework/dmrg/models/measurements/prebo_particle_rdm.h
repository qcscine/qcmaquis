/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2021 by Robin Feldmann <robinfe@phys.chem.ethz.ch>
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

#ifndef MAQUIS_DMRG_PREBO_PARTICLE_RDM_H
#define MAQUIS_DMRG_PREBO_PARTICLE_RDM_H


/* internal */
//#include "dmrg/block_matrix/symmetry/nu1.h"
#include "dmrg/models/measurement.h"
#include "dmrg/utils/checks.h"
#include "measurements_details.h"
#include "dmrg/models/generate_mpo/utils.hpp"
#include "dmrg/models/coded/models_nu1.hpp"
#include "dmrg/models/prebo/prebo_TermGenerator.hpp"
/* external */
#include <Eigen/Dense>
#include <Eigen/Eigenvalues>

namespace measurements {

    template <class Matrix, class NU1, class = void>
    class PreBOParticleRDM : public measurement<Matrix, NU1> {
        typedef typename generate_mpo::OperatorTagTerm<Matrix, NU1>::tag_type tag_type;
        typedef std::vector<tag_type> tag_vec;
    public:
        PreBOParticleRDM(BaseParameters& parms, Lattice const& lat_, const tag_vec& identities, const tag_vec& fillings,
                         std::shared_ptr<prebo::TermGenerator<Matrix, NU1>> ptr_term_generator_) : parms(parms),
                         lat(lat_), identities(identities), fillings(fillings), ptr_term_generator(ptr_term_generator_) {

        }

        ~PreBOParticleRDM()=default;

        void evaluate(MPS<Matrix, NU1> const& ket_mps, boost::optional<reduced_mps<Matrix, NU1> const&> rmps = boost::none) {
            auto vec_orbitals = lat.template get_prop<std::vector<int>>("vec_orbitals");
            auto isFermion = lat.template get_prop<std::vector<bool>>("isFermion");
            auto inv_order = lat.template get_prop<std::vector<Lattice::pos_t>>("inv_order");
            for (int nt=0; nt<vec_orbitals.size(); ++nt) {
                size_t dim = vec_orbitals.at(nt);
                Eigen::MatrixXd rdm = Eigen::MatrixXd::Zero(vec_orbitals.at(nt), vec_orbitals.at(nt));
                for (int mu = 0; mu < dim; mu++) {
                    for (int nu = 0; nu <= mu; nu++) {
                        auto terms = ptr_term_generator->generate_terms1RDM(nt, mu, nu);
                        generate_mpo::TaggedMPOMaker<Matrix, NU1> mpom(lat, identities, fillings, ptr_term_generator->getTagHandler(), terms);
                        MPO<Matrix, NU1> mpo = mpom.create_mpo();
                        MPO<Matrix, NU1> mpoc = mpo;
                        rdm(mu, nu) = maquis::real(expval(ket_mps, mpoc));
                        if (nu != mu)
                            rdm(nu, mu) = rdm(mu, nu);
                    }
                }
                std::cout << std::endl;
                std::cout << "The one-body RDM is" << std::endl;
                std::cout << rdm << std::endl;
                std::cout << std::endl;
                Eigen::SelfAdjointEigenSolver<Eigen::MatrixXd> es(rdm);
                std::cout << "The eigenvalues of the one-body RDM are:" << std::endl << es.eigenvalues() << std::endl;
                std::cout << "The trace is:      " << rdm.trace() << std::endl;
                std::cout << "The norm is:    " << norm(ket_mps) << std::endl;
                if (parms.is_set("rdmfile")) {
                    std::cout << "Writing 1-body rdm on disk ... " << std::endl;
                    std::ofstream file(parms["rdmfile"].c_str());
                    const static Eigen::IOFormat CSVFormat(Eigen::FullPrecision, Eigen::DontAlignCols, ",", "\n");
                    file << rdm.format(CSVFormat);
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
        std::shared_ptr<prebo::TermGenerator<Matrix, NU1>> ptr_term_generator;
    };
}





#endif //MAQUIS_DMRG_PREBO_PARTICLE_RDM_H
