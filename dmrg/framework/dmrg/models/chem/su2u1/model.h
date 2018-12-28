/*****************************************************************************
 *
 * QCMaquis DMRG Project
 *
 * Copyright (C) 2013 Laboratory for Physical Chemistry, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2018 by Leon Freitag <lefreita@ethz.ch>
 *
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

#ifndef QC_SU2_H
#define QC_SU2_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>
#include <boost/filesystem.hpp>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/util.h"
#include "dmrg/models/chem/parse_integrals.h"
#include "dmrg/models/chem/pg_util.h"
#include "dmrg/models/chem/su2u1/chem_helper.h"
#include "dmrg/models/chem/su2u1/term_maker.h"

template<class Matrix, class SymmGroup>
class qc_su2 : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;

    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:

    qc_su2(Lattice const & lat_, BaseParameters & parms_);

    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }

    void create_terms();

    // For this model: site_type == point group irrep
    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        return phys_indices[type];
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident[type];
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill[type];
    }

    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return chem_detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "ident_full")
            return ident_full[type];

        else if (name == "create")
            return create[type];
        else if (name == "create_fill")
            return create_fill[type];
        else if (name == "create_couple_up")
            return create_couple_up[type];
        else if (name == "create_fill_couple_down")
            return create_fill_couple_down[type];

        else if (name == "destroy")
            return destroy[type];
        else if (name == "destroy_fill")
            return destroy_fill[type];
        else if (name == "destroy_couple_up")
            return destroy_couple_up[type];
        else if (name == "destroy_fill_couple_down")
            return destroy_fill_couple_down[type];

        else if (name == "count")
            return count[type];
        else if (name == "count_fill")
            return count_fill[type];

        else if (name == "docc")
            return docc[type];
        else
            throw std::runtime_error("Operator not valid for this model: " + name);
        return 0;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }

    measurements_type measurements () const
    {
        typedef boost::tokenizer<boost::char_separator<char> > tokenizer;

        measurements_type meas;

        typedef std::vector<tag_type> tag_vec;
        typedef std::vector<tag_vec> bond_element;

        // Regexps for RDM and TDM measurement
        boost::regex expression_oneptdm("^MEASURE\\[1rdm\\]");
        boost::regex expression_twoptdm("^MEASURE\\[2rdm\\]");

        boost::regex expression_transition_oneptdm("^MEASURE\\[trans1rdm\\]");
        boost::regex expression_transition_twoptdm("^MEASURE\\[trans2rdm\\]");

        // Regexps for RDM derivative measurement
        boost::regex expression_onerdm_derivativeL("^MEASURE\\[1rdm-derivativeL\\]");
        boost::regex expression_twordm_derivativeL("^MEASURE\\[2rdm-derivativeL\\]");
        boost::regex expression_onerdm_derivativeR("^MEASURE\\[1rdm-derivativeR\\]");
        boost::regex expression_twordm_derivativeR("^MEASURE\\[2rdm-derivativeR\\]");

        // If 'Xrdm-derivative' without a suffix is specified, both left and right derivatives
        // are calculated
        boost::regex expression_onerdm_derivative_both("^MEASURE\\[1rdm-derivative\\]");
        boost::regex expression_twordm_derivative_both("^MEASURE\\[2rdm-derivative\\]");

        // Regexp for local Hamiltonian matrix elements
        boost::regex expression_local_hamiltonian("^MEASURE\\[local-hamiltonian\\]");
        // Regexp for diagonal local Hamiltonian matrix elements
        boost::regex expression_local_hamiltonian_diag("^MEASURE\\[local-hamiltonian-diag\\]");
        // or sigma vector (contracted local Hamiltonian with one MPS)
        boost::regex expression_sigma_vector("^MEASURE\\[sigma-vector\\]");

        // Regexp for Lagrange RDM update (MPS contribution to the Lagrange effective RDM in gradient calculations)
        // as for the RDM derivatives, the name without a L or R suffix means both left and right expectation values
        // are calculated
        boost::regex expression_onerdm_lagrangeL("^MEASURE\\[1rdm-lagrangeL\\]");
        boost::regex expression_twordm_lagrangeL("^MEASURE\\[2rdm-lagrangeL\\]");
        boost::regex expression_onerdm_lagrangeR("^MEASURE\\[1rdm-lagrangeR\\]");
        boost::regex expression_twordm_lagrangeR("^MEASURE\\[2rdm-lagrangeR\\]");
        boost::regex expression_onerdm_lagrange_both("^MEASURE\\[1rdm-lagrange\\]");
        boost::regex expression_twordm_lagrange_both("^MEASURE\\[2rdm-lagrange\\]");

        // Regexp for dumping a TwoSiteTensor at site X where X is specified as "MEASURE[dump-tst] = X"
        boost::regex expression_dump_tst("^MEASURE\\[dump-tst\\]");

        boost::smatch what;

        for (alps::Parameters::const_iterator it=parms.begin();it != parms.end();++it) {
            std::string lhs = it->key();

            std::string name, nameR;
            std::string bra_ckp("");
            bool expr_rdm = false, expr_rdm_derivative = false, expr_rdm_deriv_both = false,
                 expr_rdm_lagrange = false, expr_rdm_lagrange_both = false;
            std::vector<pos_t> positions;
            // Measure 1-RDM, 2-RDM, 1-TDM or 2-TDM
            // for TDMs the measurement is <bra_ckp| (operators) | this>

            if (boost::regex_match(lhs, what, expression_oneptdm)) {

                name = "oneptdm";
                expr_rdm = true;
            }

            if (boost::regex_match(lhs, what, expression_twoptdm)) {

                name = "twoptdm";
                expr_rdm = true;

            }

            if (boost::regex_match(lhs, what, expression_transition_oneptdm)) {

                name = "transition_oneptdm";
                bra_ckp = it->value();
                expr_rdm = true;
            }

            if (boost::regex_match(lhs, what, expression_transition_twoptdm)) {

                name = "transition_twoptdm";
                bra_ckp = it->value();
                expr_rdm = true;
            }

            // Measure RDM derivatives
            if (boost::regex_match(lhs, what, expression_onerdm_derivativeL)) {

                name = "onerdmderivL";
                expr_rdm_derivative = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_derivativeL)) {

                name = "twordmderivL";
                expr_rdm_derivative = true;
            }

            if (boost::regex_match(lhs, what, expression_onerdm_derivativeR)) {

                name = "onerdmderivR";
                expr_rdm_derivative = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_derivativeR)) {

                name = "twordmderivR";
                expr_rdm_derivative = true;
            }

            if (boost::regex_match(lhs, what, expression_onerdm_derivative_both)) {
                name = "onerdmderivL"; nameR = "onerdmderivR";
                expr_rdm_deriv_both = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_derivative_both)) {
                name = "twordmderivL"; nameR = "twordmderivR";
                expr_rdm_deriv_both = true;
            }

            // Measure MPS contributions to the effective RDM from Lagrange multipliers in linear response equations
            if (boost::regex_match(lhs, what, expression_onerdm_lagrangeL)) {

                name = "onerdmlagrangeL";
                expr_rdm_lagrange = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_lagrangeL)) {

                name = "twordmlagrangeL";
                expr_rdm_lagrange = true;
            }

            if (boost::regex_match(lhs, what, expression_onerdm_lagrangeR)) {

                name = "onerdmlagrangeR";
                expr_rdm_lagrange = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_lagrangeR)) {

                name = "twordmlagrangeR";
                expr_rdm_lagrange = true;
            }

            if (boost::regex_match(lhs, what, expression_onerdm_lagrange_both)) {
                name = "onerdmlagrangeL"; nameR = "onerdmlagrangeR";
                expr_rdm_lagrange_both = true;
            }

            if (boost::regex_match(lhs, what, expression_twordm_lagrange_both)) {
                name = "twordmlagrangeL"; nameR = "twordmlagrangeR";
                expr_rdm_lagrange_both = true;
            }

            if (expr_rdm)
                meas.push_back( new measurements::TaggedNRankRDM<Matrix, SymmGroup>(
                                name, lat, tag_handler, op_collection, positions, bra_ckp));

            if (expr_rdm_derivative)
                meas.push_back( new measurements::NRDMDerivative<Matrix, SymmGroup>(
                                parms["lrparam_site"], // Site for LR parameters
                                symm_traits::HasSU2<SymmGroup>(), // specialization of the constructor for SU2U1
                                name, lat, tag_handler, op_collection, positions));

            if (expr_rdm_deriv_both) {
                meas.push_back( new measurements::NRDMDerivative<Matrix, SymmGroup>(
                                parms["lrparam_site"],
                                symm_traits::HasSU2<SymmGroup>(),
                                name, lat, tag_handler, op_collection, positions));
                meas.push_back( new measurements::NRDMDerivative<Matrix, SymmGroup>(
                                parms["lrparam_site"],
                                symm_traits::HasSU2<SymmGroup>(),
                                nameR, lat, tag_handler, op_collection, positions));
            }

            if (expr_rdm_lagrange)
                meas.push_back( new measurements::NRDMLRLagrange<Matrix, SymmGroup>(
                                parms["lrparam_site"],
                                symm_traits::HasSU2<SymmGroup>(), // specialization of the constructor for SU2U1
                                it->value(), name, lat, tag_handler, op_collection, positions));

            if (expr_rdm_lagrange_both) {
                meas.push_back( new measurements::NRDMLRLagrange<Matrix, SymmGroup>(
                                parms["lrparam_site"],
                                symm_traits::HasSU2<SymmGroup>(),
                                it->value(), name, lat, tag_handler, op_collection, positions));
                meas.push_back( new measurements::NRDMLRLagrange<Matrix, SymmGroup>(
                                parms["lrparam_site"],
                                symm_traits::HasSU2<SymmGroup>(),
                                it->value(), nameR, lat, tag_handler, op_collection, positions));
            }

            // Measure Local Hamiltonian matrix elements
            if (boost::regex_match(lhs, what, expression_local_hamiltonian)) {
                name = "local_hamiltonian";
                meas.push_back( new measurements::LocalHamiltonian<Matrix, SymmGroup>(
                                name, lat, parms));
            }
            // diagonal Local Hamiltonian matrix elements
            if (boost::regex_match(lhs, what, expression_local_hamiltonian_diag)) {
                name = "local_hamiltonian_diag";
                meas.push_back( new measurements::LocalHamiltonian<Matrix, SymmGroup>(
                                name, lat, parms));
            }
            // or sigma vector
            // If we provide a valid file name (for the auxiliary MPSTensor contents) as a parameter, pass it to the constructor
            // If the file does not exist or cannot be opened, silently ignore it (TODO: Maybe this behaviour is not clean!)
            if (boost::regex_match(lhs, what, expression_sigma_vector)) {
                name = "sigma_vector";
                std::string ext_filename = boost::filesystem::exists(it->value()) ? it->value() : "";
                meas.push_back( new measurements::LocalHamiltonian<Matrix, SymmGroup>(
                                name, lat, parms, ext_filename));
            }

            // Dump the two-site tensor at site X where X is specified in lrparam_site
            if (boost::regex_match(lhs, what, expression_dump_tst)) {
                name = "tstdump";
                meas.push_back( new measurements::DumpTST<Matrix, SymmGroup>(
                                name, parms["lrparam_site"]));
            }
        }

        return meas;
    }

private:
    Lattice const & lat;
    BaseParameters & parms;
    std::vector<Index<SymmGroup> > phys_indices;

    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    std::vector<tag_type> create_fill, create, destroy_fill, destroy,
                          create_fill_couple_down, destroy_fill_couple_down,
                          create_couple_up, destroy_couple_up,
                          create_fill_count, create_count, destroy_fill_count, destroy_count,
                          count, docc, e2d, d2e, flip_S0, flip_to_S2, flip_to_S0,
                          ident, ident_full, fill, count_fill;

    typename TermMakerSU2<Matrix, SymmGroup>::OperatorCollection op_collection;

    typename SymmGroup::subcharge max_irrep;

    std::vector<op_t> generate_site_specific_ops(op_t const & op) const
    {
        PGDecorator<SymmGroup> set_symm;
        std::vector<op_t> ret;
        for (typename SymmGroup::subcharge sc=0; sc < max_irrep+1; ++sc) {
            op_t mod(set_symm(op.basis(), sc));
            mod.spin() = op.spin();
            for (std::size_t b = 0; b < op.n_blocks(); ++b)
                mod[b] = op[b];

            ret.push_back(mod);
        }
        return ret;
    }

    std::vector<tag_type> register_site_specific(std::vector<op_t> const & ops, tag_detail::operator_kind kind)
    {
        std::vector<tag_type> ret;
        for (typename SymmGroup::subcharge sc=0; sc < max_irrep+1; ++sc) {
            std::pair<tag_type, value_type> newtag = tag_handler->checked_register(ops[sc], kind);
            assert( newtag.first < tag_handler->size() );
            assert( std::abs(newtag.second - value_type(1.)) == value_type() );
            ret.push_back(newtag.first);
        }

        return ret;
    }

};

#include "dmrg/models/chem/su2u1/model.hpp"

#endif
