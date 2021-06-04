/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2021 Institute for Theoretical Physics, ETH Zurich
 *               2021- by Michele Dolfi <robinfe@ethz.ch>
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

#ifndef MAQUIS_DMRG_NU1_SYMBOLICJORDANWIGNER_HPP
#define MAQUIS_DMRG_NU1_SYMBOLICJORDANWIGNER_HPP

#ifdef HAVE_NU1

#include <dmrg/models/lattice.h>

enum class Type { Fermion, Boson };
enum class Spin { Up, Down, Zero, None };
enum class OpType { Create, Annihilate, Filling, Ident };

/**
 * @class SymbolicOperator
 * 
 * This small class represents a generic second-quantization operator.
 * This is needed to generate the proper terms based on the Jordan--Wigner transformation.
 */
class SymbolicOperator {
public:
    // Types definition
    using part_type = typename Lattice::part_type;
    using pos_t =  typename Lattice::pos_t;
private:
    // Class members
    part_type m_part_type;
    Spin m_spin;
    pos_t m_site;
    OpType m_opType;

public:

    /** @brief Constructor */
    SymbolicOperator(pos_t site, OpType opType, part_type pt, Spin spin = Spin::None)
            : m_site(site), m_opType(opType), m_part_type(pt), m_spin(spin) {}

    /**
     * @brief Copy constructor + op type
     * Note that in this case the operator type is given externally, and is not read
     * from the SymOp object.
     */
    SymbolicOperator(const SymbolicOperator &SymOp, OpType opType) {
        m_part_type = SymOp.getPartType();
        if (opType == OpType::Filling || opType == OpType::Ident)
            m_spin = Spin::None;
        else
            m_spin = SymOp.getSpin();
        m_site = SymOp.getSite();
        m_opType = opType;
    }

    /** @brief Copy constructor */
    SymbolicOperator(const SymbolicOperator &SymOp) {
        m_part_type = SymOp.getPartType();
        m_spin = SymOp.getSpin();
        m_site = SymOp.getSite();
        m_opType = SymOp.getOpType();
    }
    
    /** @brief Spin getter */
    Spin getSpin() const { return m_spin; }

    /** @brief Site getter */
    pos_t getSite() const { return m_site; }

    /** @brief Op getter */
    OpType getOpType() const { return m_opType; }

    /** @brief Particle type getter */
    part_type getPartType() const { return m_part_type; }

    /** @brief This method prints the operator in an algebraic notation */
    void print() const {
        std::string op_str;
        if (m_opType == OpType::Filling)
            op_str.append("F");
        else if (m_opType == OpType::Ident)
            op_str.append("I");
        else if (m_spin == Spin::Up || m_spin == Spin::Down)
            op_str.append("a");
        else if (m_spin == Spin::Zero)
            op_str.append("b");
        if (m_opType == OpType::Create)
            op_str.append("+");
        op_str.append("^");
        std::stringstream ss;
        ss << m_part_type;
        op_str.append(ss.str());
        op_str.append("_");
        ss.str(std::string());
        ss << m_site;
        op_str.append(ss.str());
        if (m_spin == Spin::Up)
            op_str.append("-up");
        else if (m_spin == Spin::Down)
            op_str.append("-down");
        std::cout << op_str;
    }
};


/**
 * @class SymbolicJordanWigner nu1_SymbolicJordanWigner.hpp
 * @brief This class handels the Jordan-Wigner transformation, based on the `SymbolicOperator`.
 * @note This is the most important and delicate part of the pre-BO model.
 * 
 * For a description of the algorithm, see:
 * Correlation effects in Multicomponent Quantum Chemistry (2020), Robin Feldmann, Markus Reiher group ETH Zurich, p 27
 */
class SymbolicJordanWigner {
    /** Types definition */
    using part_type = typename Lattice::part_type;
    using pos_t = typename Lattice::pos_t;
private:
    /** Class members */
    std::vector<SymbolicOperator> m_SymOpStr;
    std::vector<std::pair<part_type, pos_t>> m_nb_term;
    bool verbose = false;

public:
    /** @brief Default class constructor */
    SymbolicJordanWigner() = default;

    /** 
     * @brief Constructor from a vector of symbolic operators
     * At construction, the object stores 
     */
    explicit SymbolicJordanWigner(std::vector<SymbolicOperator> SymOpStr) {

        if (SymOpStr.size() == 0)
            throw std::runtime_error("Empty operator string in JW.");
        std::vector<part_type> vec_pt_check;
        // Identity check:
        // If the operator string contains only identities, no JW has to be applied.
        bool allIdent=true;
        for (auto const& Op : SymOpStr) {
            vec_pt_check.push_back(Op.getPartType());
            if (allIdent && Op.getOpType()!=OpType::Ident)
                allIdent=false;
        }
        if (!(std::equal(vec_pt_check.begin() + 1, vec_pt_check.end(), vec_pt_check.begin())))
            throw std::runtime_error("JW must only be applied to operator strings containing equal particle types.");

        m_SymOpStr = SymOpStr;

        if (!allIdent) {
            // Identity removal
            for (auto it=m_SymOpStr.begin(); it!=m_SymOpStr.end();) {
                if (it->getOpType() == OpType::Ident)
                    it = m_SymOpStr.erase(it);
                else
                    ++it;
            }
            for (auto const& iter : m_SymOpStr)
                m_nb_term.push_back(std::make_pair(iter.getPartType(), iter.getSite()));

            std::vector<pos_t> lindex;
            std::vector<pos_t> rindex;

            JW_generate_lr_index(m_nb_term, lindex, rindex);
            insertFillings(lindex, rindex);
        }
        if (verbose) {
            for (auto const& symOp : m_SymOpStr) {
                symOp.print();
                std::cout << " ";
            }
            std::cout << std::endl;
        }
    }

    /** @brief Getter for the symbolic operator string */
    const std::vector<SymbolicOperator> &getSymOpStr() const { return m_SymOpStr; }

private:
    /**
     * @brief Jordan--Wigner Filling
     * 
     * This method is the HEART and BRAIN of the Jordan--Wigner transformation.
     * It evaluates which operator gets assigned a filling operator on the left
     * and right side.
     */
    static void
    JW_generate_lr_index(const std::vector<std::pair<part_type, pos_t>> &nb_term,
                         std::vector<pos_t> &lindex, std::vector<pos_t> &rindex) {
        unsigned int length = nb_term.size();
        for (unsigned int i = 0; i < length; i++) {
            lindex.push_back(0);
            rindex.push_back(0);
        }
        //
        // VERY IMPORTANT:  This std::vector of maps takes care of the Jordan--Wigner
        // transformation!!
        //                  --> the procedure is very delicate and error-prone,
        //                  --> any modifications must be debugged very carefully
        //
        // This std::vector of maps stores the filling string for the relevant sites for
        // each operator. key: index, value: true if operator at position in std::vector
        // provides unassigned filling for site ´´índex´´.
        std::vector<std::map<pos_t, bool>> vec_filling_map;
        for (unsigned int i = 0; i < length; i++) {
            std::map<pos_t, bool> temp_map;
            // For each operator: loop over other operators and if index is larger,
            // then set to true. Can only be set to true once.
            for (unsigned int j = 0; j < length; j++) {
                if (temp_map.count(nb_term[j].second) == 0) {
                    if (nb_term[i].second > nb_term[j].second) {
                        temp_map.insert(std::make_pair(nb_term[j].second, true));
                    } else
                        temp_map.insert(std::make_pair(nb_term[j].second, false));
                }
            }
            vec_filling_map.push_back(temp_map);
        }

        //
        // Each operator with a larger index can only provide ONE filling operator
        // for an operator with a smaller index. This is taken care of in the
        // following. If the operator had a "true" for another operator it will be
        // set "false" if the count was raised.
        //

        // Loop from left to right in order to fill lindex.
        for (unsigned int i = 1; i < length; i++) {
            unsigned int lcount = 0;
            for (unsigned int j = 0; j < i; j++) {
                if (nb_term[i].second < nb_term[j].second &&
                    nb_term[i].first == nb_term[j].first &&
                    vec_filling_map.at(j)[nb_term[i].second]) {
                    lcount++;
                    vec_filling_map.at(j)[nb_term[i].second] = false;
                }
            }
            lindex.at(i) = lcount;
        }
        // Loop from right to left in order to fill rindex
        for (int i = length - 2; i >= 0; i--) {
            unsigned int rcount = 0;
            for (int j = length - 1; j > i; j--) {
                if (nb_term[i].second < nb_term[j].second &&
                    nb_term[i].first == nb_term[j].first &&
                    vec_filling_map.at(j)[nb_term[i].second] == true) {
                    rcount++;
                    vec_filling_map.at(j)[nb_term[i].second] = false;
                }
            }
            rindex.at(i) = rcount;
        }
        // Finally, two filling operators cancel each other out.
        for (unsigned int i = 0; i < length; i++) {
            lindex.at(i) = lindex.at(i) % 2;
            rindex.at(i) = rindex.at(i) % 2;
        }
    }

    void insertFillings(const std::vector<pos_t> &lindex,
                        const std::vector<pos_t> &rindex) {
        std::vector<SymbolicOperator> new_SymOpStr;
        size_t count = 0;
        for (auto const &symOp : m_SymOpStr) {
            // fermion
            if (symOp.getSpin() == Spin::Down ||
                symOp.getSpin() == Spin::Up) {

                // c-(d) := F a(d)
                if (symOp.getSpin() == Spin::Down &&
                    symOp.getOpType() == OpType::Annihilate)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, OpType::Filling));

                if (lindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, OpType::Filling));

                new_SymOpStr.push_back(SymbolicOperator(symOp));

                // c+(d) := a+(d) F
                if (symOp.getSpin() == Spin::Down &&
                    symOp.getOpType() == OpType::Create)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, OpType::Filling));

                if (rindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, OpType::Filling));
            }
            // boson
            else {
                new_SymOpStr.push_back(SymbolicOperator(symOp));
            }
            count++;
        }

        // Cleaning of the JW string.
        // Probs go to Alberto :)
        for (auto i = new_SymOpStr.begin(); i != new_SymOpStr.end();) {
            auto n = std::next(i);
            if (n == new_SymOpStr.end())
                break;
            if (i->getPartType() == n->getPartType() &&
                i->getSite() == n->getSite() && i->getOpType() == n->getOpType() &&
                i->getOpType() == OpType::Filling) {
                i = new_SymOpStr.erase(i);
                i = new_SymOpStr.erase(i);
            } else {
                i++;
            }
        }
        m_SymOpStr = new_SymOpStr;
    }
};

#endif

#endif //MAQUIS_DMRG_NU1_SYMBOLICJORDANWIGNER_HPP