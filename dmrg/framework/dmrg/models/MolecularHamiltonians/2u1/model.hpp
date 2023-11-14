/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2021- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef QC_MODEL_HPP
#define QC_MODEL_HPP

#include "dmrg/models/JordanWignerManager.h"
#include "normal_order.h"

using Hamiltonian = chem::Hamiltonian;
using HamiltonianTransformation = chem::HamiltonianTransformation;

template <class Matrix, class SymmGroup, Hamiltonian HamiltonianType, HamiltonianTransformation Transcorrelated>
qc_model<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::qc_model(Lattice const & lat_, BaseParameters & parms_)
    : lat(lat_), parms(parms_), tag_handler(new table_type()), isQuantumComputingFormat(false)
{
    // Types definition
    typedef typename SymmGroup::subcharge subcharge;

    // Parameter parsing
    if (isTranscorrelated_ && parms.is_set("transcorrelated_quantum_computing_format")) {
        if (parms["transcorrelated_quantum_computing_format"] == "yes") {
            maquis::cout << " Activating transcorrelated quantum computing format" << std::endl;
            isQuantumComputingFormat = true;
        }
    }

    if (!isTranscorrelated_ && parms.is_set("quantum_computing_format")) {
        if (parms["quantum_computing_format"] == "yes") {
            maquis::cout << " Activating conventional quantum computing format" << std::endl;
            isQuantumComputingFormat = true;
        }
    }

    // find the highest irreducible representation number
    // used to generate ops for all irreps 0..max_irrep
    max_irrep = 0;
    for (pos_t p = 0; p < lat.size(); ++p)
        max_irrep = (lat.get_prop<typename SymmGroup::subcharge>("type", p) > max_irrep)
                    ? lat.get_prop<typename SymmGroup::subcharge>("type", p) : max_irrep;

    typename SymmGroup::charge A(0), B(0), C(0), D(1);
    B[0] = 1;
    C[1] = 1;

    for (subcharge irr = 0; irr <= max_irrep; ++irr) {
        Index<SymmGroup> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(B, irr), 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(C, irr), 1));
        phys.insert(std::make_pair(D, 1));

        phys_indices.push_back(phys);
    }

    op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
            count_up_op, count_down_op, count_up_down_op, docc_op, e2d_op, d2e_op,
            d2u_op, u2d_op, create_down_for_meas_op, destroy_down_for_meas_op,
            ident_op, fill_op;

    ident_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_op.insert_block(Matrix(1, 1, 1), B, B);
    ident_op.insert_block(Matrix(1, 1, 1), C, C);
    ident_op.insert_block(Matrix(1, 1, 1), D, D);

    create_up_op.insert_block(Matrix(1, 1, 1), A, B);
    create_up_op.insert_block(Matrix(1, 1, 1), C, D);
    create_down_op.insert_block(Matrix(1, 1, 1), A, C);
    create_down_op.insert_block(Matrix(1, 1, 1), B, D);

    destroy_up_op.insert_block(Matrix(1, 1, 1), B, A);
    destroy_up_op.insert_block(Matrix(1, 1, 1), D, C);
    destroy_down_op.insert_block(Matrix(1, 1, 1), C, A);
    destroy_down_op.insert_block(Matrix(1, 1, 1), D, B);

    count_up_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_op.insert_block(Matrix(1, 1, 1), D, D);

    count_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_down_op.insert_block(Matrix(1, 1, 1), D, D);

    count_up_down_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_up_down_op.insert_block(Matrix(1, 1, 2), D, D);

    docc_op.insert_block(Matrix(1, 1, 1), D, D);

    e2d_op.insert_block(Matrix(1, 1, 1), A, D);
    d2e_op.insert_block(Matrix(1, 1, 1), D, A);

    fill_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_op.insert_block(Matrix(1, 1, -1), B, B);
    fill_op.insert_block(Matrix(1, 1, -1), C, C);
    fill_op.insert_block(Matrix(1, 1, 1), D, D);

    op_t tmp;

    // TODO ALB FOR NOW KEPT, BUT THIS SHOULD GO!!
    gemm(fill_op, create_down_op, tmp);
    create_down_for_meas_op = tmp;
    gemm(destroy_down_op, fill_op, tmp);
    destroy_down_for_meas_op = tmp;

    /// stknecht: needed for special 1-TDMs
    gemm(destroy_down_op, create_up_op, d2u_op); // S_plus
    gemm(destroy_up_op, create_down_op, u2d_op); // S_minus

    // only effective if point group symmetry is active, need to adapt operators to different irreps
#define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);
    GENERATE_SITE_SPECIFIC(ident_op)
    GENERATE_SITE_SPECIFIC(fill_op)
    GENERATE_SITE_SPECIFIC(create_up_op)
    GENERATE_SITE_SPECIFIC(create_down_op)
    GENERATE_SITE_SPECIFIC(create_down_for_meas_op)
    GENERATE_SITE_SPECIFIC(destroy_up_op)
    GENERATE_SITE_SPECIFIC(destroy_down_op)
    GENERATE_SITE_SPECIFIC(destroy_down_for_meas_op)
    GENERATE_SITE_SPECIFIC(count_up_op)
    GENERATE_SITE_SPECIFIC(count_down_op)
    GENERATE_SITE_SPECIFIC(e2d_op)
    GENERATE_SITE_SPECIFIC(d2e_op)
    GENERATE_SITE_SPECIFIC(docc_op)
    GENERATE_SITE_SPECIFIC(count_up_down_op)
    GENERATE_SITE_SPECIFIC(d2u_op)
    GENERATE_SITE_SPECIFIC(u2d_op)
#undef GENERATE_SITE_SPECIFIC

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

#define REGISTER(op, kind) op = this->register_site_specific(op ## _ops, kind);
    REGISTER(ident, tag_detail::bosonic)
    REGISTER(fill, tag_detail::bosonic)
    REGISTER(create_up, tag_detail::fermionic)
    REGISTER(create_down, tag_detail::fermionic)
    REGISTER(create_down_for_meas, tag_detail::fermionic)
    REGISTER(destroy_up, tag_detail::fermionic)
    REGISTER(destroy_down, tag_detail::fermionic)
    REGISTER(destroy_down_for_meas, tag_detail::fermionic)
    REGISTER(count_up, tag_detail::bosonic)
    REGISTER(count_down, tag_detail::bosonic)
    REGISTER(e2d, tag_detail::bosonic)
    REGISTER(d2e, tag_detail::bosonic)
    REGISTER(docc, tag_detail::bosonic)
    REGISTER(count_up_down, tag_detail::bosonic)
    REGISTER(d2u, tag_detail::bosonic)
    REGISTER(u2d, tag_detail::bosonic)
#undef REGISTER

    //**********************************************************************
    std::pair<std::vector<tag_type>, std::vector<value_type> > cutf = tag_handler->get_product_tags(create_up, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > cdtf = tag_handler->get_product_tags(
            create_down_for_meas, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ftdu = tag_handler->get_product_tags(fill, destroy_up);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ftdd = tag_handler->get_product_tags(fill,
                                                                                                    destroy_down_for_meas);
    std::pair<std::vector<tag_type>, std::vector<value_type> > cund = tag_handler->get_product_tags(create_up,
                                                                                                    count_down);
    std::pair<std::vector<tag_type>, std::vector<value_type> > dund = tag_handler->get_product_tags(destroy_up,
                                                                                                    count_down);
    std::pair<std::vector<tag_type>, std::vector<value_type> > cdnu = tag_handler->get_product_tags(
            create_down_for_meas, count_up);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ddnu = tag_handler->get_product_tags(
            destroy_down_for_meas, count_up);
    std::pair<std::vector<tag_type>, std::vector<value_type> > cundtf = tag_handler->get_product_tags(cund.first, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ftdund = tag_handler->get_product_tags(fill, dund.first);
    std::pair<std::vector<tag_type>, std::vector<value_type> > cdnutf = tag_handler->get_product_tags(cdnu.first, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ftddnu = tag_handler->get_product_tags(fill, ddnu.first);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ddcu = tag_handler->get_product_tags(
            destroy_down_for_meas, create_up);
    std::pair<std::vector<tag_type>, std::vector<value_type> > ducd = tag_handler->get_product_tags(destroy_up,
                                                                                                    create_down_for_meas);

    // Note that the Hermitian pairs are registered only if the Hamiltonian is Hermitean.
    // TODO: In principle, also for the transcorrelated case the registration of the hermitean pairs should
    //       work, needs more testing to understand why it does not work.
    if (!isTranscorrelated_) {
        int numberOfTypes = create_up.size();
        for (int opType = 0; opType < numberOfTypes; opType++) {
            tag_handler->hermitian_pair(create_up[opType], destroy_up[opType]);
            tag_handler->hermitian_pair(create_down[opType], destroy_down[opType]);
            tag_handler->hermitian_pair(create_down_for_meas[opType], destroy_down_for_meas[opType]);
            tag_handler->hermitian_pair(cutf.first[opType], ftdu.first[opType]);
            tag_handler->hermitian_pair(cdtf.first[opType], ftdd.first[opType]);
            tag_handler->hermitian_pair(e2d[opType], d2e[opType]);
            tag_handler->hermitian_pair(cund.first[opType], dund.first[opType]);
            tag_handler->hermitian_pair(cdnu.first[opType], ddnu.first[opType]);
            tag_handler->hermitian_pair(cundtf.first[opType], ftdund.first[opType]);
            tag_handler->hermitian_pair(cdnutf.first[opType], ftddnu.first[opType]);
            tag_handler->hermitian_pair(ddcu.first[opType], ducd.first[opType]);
        }
    }
    if (isTranscorrelated_)
        maquis::cout << "Transcorrelated Hamiltonian modality activated" << std::endl;
}

/** @brief Create the Hamiltonian terms */
template<class Matrix, class SymmGroup, Hamiltonian HamiltonianType, HamiltonianTransformation Transcorrelated>
void qc_model<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::create_terms() {
    bool is_normal_ordered = (parms["transcorrelated_3body_normal_ordered"] == "yes");

    if(isTranscorrelated_ && is_normal_ordered) {
        create_terms_normal_ordered();
    } else {
        create_terms_not_normal_ordered();
    }
}

/** @brief Create the Hamiltonian terms */
template<class Matrix, class SymmGroup, Hamiltonian HamiltonianType, HamiltonianTransformation Transcorrelated>
void qc_model<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::create_terms_not_normal_ordered() {
    // Generates the data required to form the Hamiltonian
    auto jw = JordanWignerHandler<Matrix, SymmGroup>(lat, fill, create_up, create_down, destroy_up, destroy_down);
    MapOfOperatorsType mapOfOperators;
    bool isTcAndQuantum = isQuantumComputingFormat && isTranscorrelated_;
    chem::detail::ChemHelper<Matrix, SymmGroup, HamiltonianType, Transcorrelated> term_assistant(parms, lat, ident, fill, tag_handler, !isTcAndQuantum);
    auto& matrix_elements = term_assistant.getMatrixElements();
    // Tmp objects.
    std::vector<OperatorType> oneBodyVec1 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> oneBodyVec2 = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
    std::vector<std::vector<OperatorType> > oneBodyElementaryOperators = {oneBodyVec1, oneBodyVec2};

    std::vector<OperatorType> twoBodyVec1 = {OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                             OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
    std::vector<OperatorType> twoBodyVec2 = {OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                             OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
    std::vector<OperatorType> twoBodyVec3 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                             OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> twoBodyVec4 = {OperatorType::CreateBeta, OperatorType::CreateBeta,
                                             OperatorType::DestroyBeta, OperatorType::DestroyBeta};
    std::vector<std::vector<OperatorType> > twoBodyElementaryOperators = {twoBodyVec1, twoBodyVec2, twoBodyVec3, twoBodyVec4};

    bool normal_ordered_integral = (parms["normal_ordered_integral_file"] == "yes") && isTranscorrelated_;
    std::unordered_set<std::size_t> hole_states;
    if(normal_ordered_integral) {
        std::vector<std::size_t> hole_state_vec = parms["normal_ordered_hole_states"];
        hole_states.insert(hole_state_vec.begin(), hole_state_vec.end());
        std::cout << "Reading integral in normally ordered form, with reference states ";
        for(auto i : hole_states) {
            std::cout << i << ",";
        }
        std::cout << std::endl;
    }

    for (std::size_t iElement = 0; iElement < matrix_elements.size(); iElement++) {
        int i = term_assistant.idx(iElement, 0);
        int j = term_assistant.idx(iElement, 1);
        int k = term_assistant.idx(iElement, 2);
        int l = term_assistant.idx(iElement, 3);
        int m = -1;
        int n = -1;
        auto matrixElement = static_cast<value_type>(matrix_elements[iElement]);
        if (isTranscorrelated_) {
            m = term_assistant.idx(iElement, 4);
            n = term_assistant.idx(iElement, 5);
        }
        // Core electrons energy
        if (i == -1 && j == -1 && k == -1 && l == -1 && m == -1 && n == -1) {
            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back(std::make_pair(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]));

            this->terms_.push_back(term);

        }
        // One-body contribution
        else if (k == -1 && l == -1 && m==-1 && n==-1) {
            std::vector< std::array<int, 2> > posVector = isTranscorrelated_ ? std::vector<std::array<int, 2>>({std::array<int, 2>({i, j})})
                                                                             : TermMaker<Matrix, SymmGroup>::generateTwofoldSymmetricIndex(i, j);
            for (auto& iOp: oneBodyElementaryOperators) {
                for (auto& iTerm: posVector) {
                    std::vector< pos_t > localPosVector = { iTerm[0], iTerm[1] };
                    if(!normal_ordered_integral) {
                        auto term = jw.getTerm(localPosVector, iOp, tag_handler, true, matrixElement);
                        addTerm(mapOfOperators, term);
                    } else {
                        auto noOp(iOp);
                        double sign = applyNormalOrdering(localPosVector, noOp, hole_states);
                        auto term = jw.getTerm(localPosVector, noOp, tag_handler, true, sign * matrixElement);
                        addTerm(mapOfOperators, term);
                    }
                }
            }
        }

        // Two-body contribution
        else if (m==-1 && n==-1) {
            std::vector< OperatorType > opVector1, opVector2, opVector3, opVector4;
            if (isQuantumComputingFormat) {
                opVector1 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha, OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                opVector2 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha, OperatorType::CreateBeta, OperatorType::DestroyBeta};
                opVector3 = {OperatorType::CreateBeta, OperatorType::DestroyBeta, OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                opVector4 = {OperatorType::CreateBeta, OperatorType::DestroyBeta, OperatorType::CreateBeta, OperatorType::DestroyBeta};
            }
            else {
                opVector1 = {OperatorType::CreateAlpha, OperatorType::CreateBeta, OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
                opVector2 = {OperatorType::CreateBeta, OperatorType::CreateAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
                opVector3 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
                opVector4 = {OperatorType::CreateBeta, OperatorType::CreateBeta, OperatorType::DestroyBeta, OperatorType::DestroyBeta};
            }
            std::vector< std::vector< OperatorType > > twoBodyElementaryOperators = { opVector1, opVector2, opVector3, opVector4 };
            std::vector< std::array<int, 4> > tmp;
            if (isTranscorrelated_)
                tmp = TermMaker<Matrix, SymmGroup>::generateTwofoldSymmetricIndex(i, j, k, l);
            else
                tmp = TermMaker<Matrix, SymmGroup>::generateEightfoldSymmetricIndex(i, j, k, l);
            //
            for (auto& iOp: twoBodyElementaryOperators) {
                for (auto& iTerm: tmp) {
                    auto posVector = (isQuantumComputingFormat) ? std::vector< pos_t >{ iTerm[0], iTerm[1], iTerm[2], iTerm[3] }
                                                                : std::vector< pos_t >{ iTerm[0], iTerm[2], iTerm[3], iTerm[1] };
                    if(!normal_ordered_integral) {
                        auto term = jw.getTerm(posVector, iOp, tag_handler, true, matrixElement / 2.);
                        if (term.size() > 0) {
                          addTerm(mapOfOperators, term);
                        }
                    } else {
                        if (!(posVector[0] == posVector[1] && iOp[0] == iOp[1]) &&
                            !(posVector[2] == posVector[3] && iOp[2] == iOp[3])) {
                            auto noOp(iOp);
                            double sign = applyNormalOrdering(posVector, noOp, hole_states);
                            auto term = jw.getTerm(posVector, noOp, tag_handler, true, sign * matrixElement / 2.);
                            if (term.size() > 0) {
                              addTerm(mapOfOperators, term);
                            }
                        }
                    }

                }
            }
        }
            // Three-body contribution
        else {
            if (isTranscorrelated_) {
                // No normal ordering
                if (parms["transcorrelated_3body"] == "yes") {
                    std::set<int> nonEqualIndices{i, j, k, l, m, n};
                    int couplingDegree = nonEqualIndices.size();
                    int maxDegree = parms["transcorrelated_3body_max_coupling"];
                    if (couplingDegree <= maxDegree) {
                        std::vector< std::array<int, 6> > tmp;
                        if (isQuantumComputingFormat)
                            tmp = std::vector<std::array<int, 6>>({std::array<int, 6>({i, j, k, l, m, n})});
                        else
                            tmp = TermMaker<Matrix, SymmGroup>::generateThreeBodySymmetricIndex(i, j, k, l, m, n);
                        std::vector< OperatorType > opVector1, opVector2, opVector3, opVector4, opVector5, opVector6, opVector7, opVector8;
                        value_type scalingFactor = (isQuantumComputingFormat) ? 1. : -1./6.;
                        if (isQuantumComputingFormat) {
                            // aaa
                            opVector1 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                            // aab
                            opVector2 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta};
                            // abb
                            opVector3 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta};
                            // aba
                            opVector4 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                            // baa
                            opVector5 = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                            // bab
                            opVector6 = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta};
                            // bbb
                            opVector7 = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta};
                            // bba
                            opVector8 = {OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateBeta,  OperatorType::DestroyBeta,
                                         OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
                        }
                        else {
                            opVector1 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                         OperatorType::DestroyAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
                            opVector2 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                         OperatorType::DestroyBeta, OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
                            opVector3 = {OperatorType::CreateAlpha, OperatorType::CreateBeta, OperatorType::CreateBeta,
                                         OperatorType::DestroyBeta, OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
                            opVector4 = {OperatorType::CreateAlpha, OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                         OperatorType::DestroyAlpha, OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
                            opVector5 = {OperatorType::CreateBeta, OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                         OperatorType::DestroyAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
                            opVector6 = {OperatorType::CreateBeta, OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                         OperatorType::DestroyBeta, OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
                            opVector7 = {OperatorType::CreateBeta, OperatorType::CreateBeta, OperatorType::CreateBeta,
                                         OperatorType::DestroyBeta, OperatorType::DestroyBeta, OperatorType::DestroyBeta};
                            opVector8 = {OperatorType::CreateBeta, OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                         OperatorType::DestroyAlpha, OperatorType::DestroyBeta, OperatorType::DestroyBeta};
                        }
                        std::vector< std::vector< OperatorType > > threeBodyElementaryOperators = { opVector1, opVector2, opVector3, opVector4,
                                                                                                    opVector5, opVector6, opVector7, opVector8 };
                        for (auto& iOp: threeBodyElementaryOperators) {
                            for (auto& iTerm: tmp) {
                                auto posVector = (isQuantumComputingFormat) ? std::vector< pos_t >{ iTerm[0], iTerm[1], iTerm[2], iTerm[3], iTerm[4], iTerm[5] }
                                                                            : std::vector< pos_t >{ iTerm[0], iTerm[2], iTerm[4], iTerm[5], iTerm[3], iTerm[1] };
                                if(!normal_ordered_integral) {
                                        auto term = jw.getTerm(posVector, iOp, tag_handler, true, - matrixElement / 6.);
                                        addTerm(mapOfOperators, term);
                                } else {
                                    if (!(posVector[0] == posVector[1] && iOp[0] == iOp[1]) &&
                                        !(posVector[0] == posVector[2] && iOp[0] == iOp[2]) &&
                                        !(posVector[1] == posVector[2] && iOp[1] == iOp[2]) &&
                                        !(posVector[4] == posVector[5] && iOp[4] == iOp[5]) &&
                                        !(posVector[3] == posVector[5] && iOp[3] == iOp[5]) &&
                                        !(posVector[4] == posVector[3] && iOp[4] == iOp[3])) {
                                            auto noOp(iOp);
                                            double sign = applyNormalOrdering(posVector, noOp, hole_states);
                                            auto term = jw.getTerm(posVector, noOp, tag_handler, true, -sign * matrixElement/ 6.);
                                            addTerm(mapOfOperators, term);
                                    }
                                }
                            }
                        }
                    }
                }

            } else {
                throw std::runtime_error("Three-body term for non transcorrelated Hamiltonian not yet available");
            }
        }
    }

    for (const auto &idx: mapOfOperators)
        this->terms_.push_back(term_descriptor(idx.first, idx.second, true));
    // Registers all Hermitian conjugate
    /*
    int originalSize = tag_handler->total_size();
    for (int iTag = 0; iTag < originalSize; iTag++) {
        auto fermType = (tag_handler->is_fermionic(iTag)) ? tag_detail::fermionic : tag_detail::bosonic;
        auto resPairCC = tag_handler->checked_register(adjoint(tag_handler->get_op(iTag)), fermType);
        if (resPairCC.first >= originalSize)
            std::cout << "Registered new Hermitian Conjugate operator" << std::endl; 202
        if (iTag < resPairCC.first) {
            if (std::abs(resPairCC.second-1.) < 1.0E-16) {
                tag_handler->hermitian_pair(iTag, resPairCC.first);
            }
            else {
                auto tmp = tag_handler->register_op(adjoint(tag_handler->get_op(iTag)), fermType);
                tag_handler->hermitian_pair(iTag, tmp);
            }
        }
    }
    */
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms" << std::endl;

}

/** @brief Create the Hamiltonian terms */
template<class Matrix, class SymmGroup, Hamiltonian HamiltonianType, HamiltonianTransformation Transcorrelated>
void qc_model<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::create_terms_normal_ordered() {
    // Generates the data required to form the Hamiltonian
    auto jw = JordanWignerHandler<Matrix, SymmGroup>(lat, fill, create_up, create_down, destroy_up, destroy_down);
    MapOfOperatorsType mapOfOperators;
    chem::detail::ChemHelper<Matrix, SymmGroup, HamiltonianType, Transcorrelated> term_assistant(parms, lat, ident,
                                                                                                 fill, tag_handler);
    auto &matrix_elements = term_assistant.getMatrixElements();
    // Tmp objects.
    std::vector<OperatorType> oneBodyVec1 = {OperatorType::CreateAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> oneBodyVec2 = {OperatorType::CreateBeta, OperatorType::DestroyBeta};
    std::vector<std::vector<OperatorType> > oneBodyElementaryOperators = {oneBodyVec1, oneBodyVec2};

    std::vector<OperatorType> twoBodyVec1 = {OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                             OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
    std::vector<OperatorType> twoBodyVec2 = {OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                             OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
    std::vector<OperatorType> twoBodyVec3 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                             OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> twoBodyVec4 = {OperatorType::CreateBeta, OperatorType::CreateBeta,
                                             OperatorType::DestroyBeta, OperatorType::DestroyBeta};
    std::vector<std::vector<OperatorType> > twoBodyElementaryOperators = {twoBodyVec1, twoBodyVec2, twoBodyVec3,
                                                                          twoBodyVec4};

    std::vector<OperatorType> threeBodyVec1 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                               OperatorType::DestroyAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> threeBodyVec2 = {OperatorType::CreateAlpha, OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                               OperatorType::DestroyBeta, OperatorType::DestroyAlpha, OperatorType::DestroyAlpha};
    std::vector<OperatorType> threeBodyVec3 = {OperatorType::CreateAlpha, OperatorType::CreateBeta, OperatorType::CreateBeta,
                                               OperatorType::DestroyBeta, OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
    std::vector<OperatorType> threeBodyVec4 = {OperatorType::CreateAlpha, OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                               OperatorType::DestroyAlpha, OperatorType::DestroyBeta, OperatorType::DestroyAlpha};
    std::vector<OperatorType> threeBodyVec5 = {OperatorType::CreateBeta, OperatorType::CreateAlpha, OperatorType::CreateAlpha,
                                               OperatorType::DestroyAlpha, OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
    std::vector<OperatorType> threeBodyVec6 = {OperatorType::CreateBeta, OperatorType::CreateAlpha, OperatorType::CreateBeta,
                                               OperatorType::DestroyBeta, OperatorType::DestroyAlpha, OperatorType::DestroyBeta};
    std::vector<OperatorType> threeBodyVec7 = {OperatorType::CreateBeta, OperatorType::CreateBeta, OperatorType::CreateBeta,
                                               OperatorType::DestroyBeta, OperatorType::DestroyBeta, OperatorType::DestroyBeta};
    std::vector<OperatorType> threeBodyVec8 = {OperatorType::CreateBeta, OperatorType::CreateBeta, OperatorType::CreateAlpha,
                                               OperatorType::DestroyAlpha, OperatorType::DestroyBeta,OperatorType::DestroyBeta};
    std::vector<std::vector<OperatorType> > threeBodyElementaryOperators = {threeBodyVec1, threeBodyVec2,
                                                                            threeBodyVec3, threeBodyVec4,
                                                                            threeBodyVec5, threeBodyVec6,
                                                                            threeBodyVec7, threeBodyVec8};
    bool added_0b = false;


    // Find hole states for normal ordering
    int n_electrons = parms["nelec"];
    int n_hole_states = (n_electrons + 1) / 2;
    std::unordered_set<std::size_t> hole_states;
    std::vector<std::pair<std::size_t, typename chem::detail::ChemHelper<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::value_type>> diagonals;
    for (std::size_t iElement = 0; iElement < matrix_elements.size(); iElement++) {
        if (term_assistant.idx(iElement, 2) == -1 && term_assistant.idx(iElement, 3) == -1
            && term_assistant.idx(iElement, 4) == -1 && term_assistant.idx(iElement, 5) == -1
            && term_assistant.idx(iElement, 0) == term_assistant.idx(iElement, 1) && term_assistant.idx(iElement, 0) != -1) { // Check if diagonal
            diagonals.push_back({iElement, matrix_elements[iElement]});
        }
    }
    std::sort(diagonals.begin(), diagonals.end(),
                      [](auto a, auto b) {
                          return a.second < b.second;
                      });
    for (auto it = diagonals.begin(); it < diagonals.begin() + n_hole_states; ++it)
        hole_states.insert(term_assistant.idx(it->first, 0));

    std::cout << "Normal ordering according to hole states ";
    for(auto h : hole_states) {
        std::cout << h << ", ";
    }
    std::cout << std::endl;

    // Normal ordering helper that generates indices for normal ordering
    NormalOrderingHelper<Matrix, SymmGroup, HamiltonianType, Transcorrelated> no_helper(term_assistant);

    for (std::size_t iElement = 0; iElement < matrix_elements.size(); iElement++) {
        int i = term_assistant.idx(iElement, 0);
        int j = term_assistant.idx(iElement, 1);
        int k = term_assistant.idx(iElement, 2);
        int l = term_assistant.idx(iElement, 3);
        int m = term_assistant.idx(iElement, 4);
        int n = term_assistant.idx(iElement, 5);
        auto matrixElement = static_cast<value_type>(matrix_elements[iElement]);

        // Core electrons energy
        if (i == -1 && j == -1 && k == -1 && l == -1 && m == -1 && n == -1) {
            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back(std::make_pair(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]));

            // NO zero body contribution
            term.coeff += no_helper.getNoZeroBodyContribution(term_assistant, matrix_elements.size(), parms["L"], hole_states);

            this->terms_.push_back(term);

            added_0b = true;
        }
            // One-body contribution
        else if (k == -1 && l == -1 && m == -1 && n == -1) {
            std::vector<std::array<int, 2> > posVector = std::vector<std::array<int, 2>>({std::array<int, 2>({i, j})});
            for (auto &iOp: oneBodyElementaryOperators) {
                for (auto &iTerm: posVector) {
                    std::vector<pos_t> localPosVector = {iTerm[0], iTerm[1]};
                    auto noOp(iOp);
                    double sign = applyNormalOrdering(localPosVector, noOp, hole_states);
                    auto term = jw.getTerm(localPosVector, noOp, tag_handler, true, sign * matrixElement);
                    addTerm(mapOfOperators, term);
                }
            }
        }
            // Normal ordered 2B contribution of 2B operator
        else if (m == -1 && n == -1) {
            std::vector<std::array<int, 4> > tmp = TermMaker<Matrix, SymmGroup>::generateTwofoldSymmetricIndex(i, j, k, l);
            for (auto &iOp: twoBodyElementaryOperators) {
                for (auto &iTerm: tmp) {
                    std::vector<pos_t> posVector = {iTerm[0], iTerm[2], iTerm[3], iTerm[1]};
                    if (!(posVector[0] == posVector[1] && iOp[0] == iOp[1]) &&
                        !(posVector[2] == posVector[3] && iOp[2] == iOp[3])) {
                        auto noOp(iOp);
                        double sign = applyNormalOrdering(posVector, noOp, hole_states);
                        auto term = jw.getTerm(posVector, noOp, tag_handler, true, sign * matrixElement / 2.);
                        addTerm(mapOfOperators, term);
                    }
                }
            }
        } else {
            // Normal ordered 3B contribution of 3B operator
            if (parms["transcorrelated_3body"] == "yes") {
                std::vector<std::array<int, 6> > tmp = TermMaker<Matrix, SymmGroup>::generateThreeBodySymmetricIndex(i, j, k, l, m, n);

                for (auto &iOp: threeBodyElementaryOperators) {
                    for (auto &iTerm: tmp) {
                        std::vector<pos_t> posVector = {iTerm[0], iTerm[2], iTerm[4], iTerm[5], iTerm[3],
                                                        iTerm[1]};
                        if (!(posVector[0] == posVector[1] && iOp[0] == iOp[1]) &&
                            !(posVector[0] == posVector[2] && iOp[0] == iOp[2]) &&
                            !(posVector[1] == posVector[2] && iOp[1] == iOp[2]) &&
                            !(posVector[4] == posVector[5] && iOp[4] == iOp[5]) &&
                            !(posVector[3] == posVector[5] && iOp[3] == iOp[5]) &&
                            !(posVector[4] == posVector[3] && iOp[4] == iOp[3])) {
                            auto noOp(iOp);
                            double sign = applyNormalOrdering(posVector, noOp, hole_states);
                            auto term = jw.getTerm(posVector, noOp, tag_handler, true,
                                                   -sign * matrixElement / 6.);
                            addTerm(mapOfOperators, term);
                        }
                    }
                }
            }
        }
    }

    // Normal ordered 2B contribution
    std::unordered_map<std::tuple<int, int, int, int>, double, intTupleHash> twoBodyIndices;
    twoBodyIndices = no_helper.getNoTwoBodyCoefficients(term_assistant, matrix_elements.size(), parms["L"],
                                                        hole_states);

    int cnt = 0;
    for (auto &iOp: twoBodyElementaryOperators) {
        for (auto p: twoBodyIndices) {
            std::vector<int> posVector = {std::get<0>(p.first), std::get<2>(p.first), std::get<3>(p.first),
                                          std::get<1>(p.first)};
            if (!(posVector[0] == posVector[1] && iOp[0] == iOp[1]) &&
                !(posVector[2] == posVector[3] && iOp[2] == iOp[3])) {
                if(cnt < 20)
                  std::cout << "Normal ordered " << posVector[0] << " " << posVector[1] <<" " << posVector[2] << " " << posVector[3] << ": " << std::setprecision(16) << p.second << std::endl;
                ++cnt;
                auto noOp(iOp);
                double sign = applyNormalOrdering(posVector, noOp, hole_states);
                auto term = jw.getTerm(posVector, noOp, tag_handler, true, sign * p.second);
                addTerm(mapOfOperators, term);
            }
        }
    }

    // Normal ordered 1B contribution
    cnt = 0;
    std::unordered_map<std::pair<int, int>, double, intPairHash> oneBodyIndices;
    oneBodyIndices = no_helper.getNoOneBodyCoefficients(term_assistant, matrix_elements.size(), parms["L"], hole_states);
    for (auto &iOp: oneBodyElementaryOperators) {
        for (auto p: oneBodyIndices) {
            std::vector<int> posVector = {p.first.first, p.first.second};
            auto noOp(iOp);
            if(cnt < 20)
                std::cout << "Normal ordered " << posVector[0] << " " << posVector[1] << ": " << std::setprecision(16) <<  p.second << std::endl;
            ++cnt;
            double sign = applyNormalOrdering(posVector, noOp, hole_states);
            auto term = jw.getTerm(posVector, noOp, tag_handler, true, sign * p.second);
            //std::cout << term << '\t';
            addTerm(mapOfOperators, term);
        }
    }

    // Normal ordered 0B contribution
    if (!added_0b) {
        term_descriptor term;
        term.push_back(std::make_pair(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]));
        term.coeff = no_helper.getNoZeroBodyContribution(term_assistant, matrix_elements.size(), parms["L"], hole_states);
        std::cout << "Adding 0B " << std::setprecision(16) << term.coeff << std::endl;
        this->terms_.push_back(term);
    }

    for (const auto &idx: mapOfOperators)
        this->terms_.push_back(term_descriptor(idx.first, idx.second, true));

    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms" << std::endl;
}

/** @brief Adds an operator to the underyling operator map */
template<class Matrix, class SymmGroup, Hamiltonian HamiltonianType, HamiltonianTransformation Transcorrelated>
void qc_model<Matrix, SymmGroup, HamiltonianType, Transcorrelated>::addTerm(MapOfOperatorsType &mapOfOperators,
                                                                            const term_descriptor &term) const {
    //
    if (term.size() != 0) {
        if (mapOfOperators.find(term.getBase()) == mapOfOperators.end()) {
            mapOfOperators.insert({term.getBase(), term.coeff });
        }
        else {
            mapOfOperators[term.getBase()] += term.coeff;
        }
    }
}

#endif
