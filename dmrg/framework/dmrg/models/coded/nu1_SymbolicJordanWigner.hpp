//
// Created by robin on 11.05.21.
//

#ifndef MAQUIS_DMRG_NU1_SYMBOLICJORDANWIGNER_HPP
#define MAQUIS_DMRG_NU1_SYMBOLICJORDANWIGNER_HPP

#if defined(HAVE_NU1)
#include <dmrg/models/lattice.h>

/**
 * @class SymbolicOperator nu1_nBodyTerm.hpp
 * @brief Symbolic algebra implementation of a 2nd quantization operator. Needed for Jordan--Wigner transformation.
 */
class SymbolicOperator {
    // Types definition
    typedef typename Lattice::part_type
            part_type; // unsigned integer denoting the particle type
    typedef typename Lattice::pos_t pos_t;

public:
    enum m_SPIN { Up, Down, Zero, None };
    enum m_OPTYPES { Create, Annihilate, Filling, Ident };

private:
    // +---------+
    // | Members |
    // +---------+
    part_type m_part_type;
    m_SPIN m_spin;
    pos_t m_site;
    m_OPTYPES m_opType;

public:
    // Constructor
    SymbolicOperator(pos_t site, m_OPTYPES opType, part_type pt,
                     m_SPIN spin = None)
            : m_site(site), m_opType(opType), m_part_type(pt), m_spin(spin) {}

    // Constructor
    SymbolicOperator(const SymbolicOperator &SymOp, m_OPTYPES opType) {
        m_part_type = SymOp.getPartType();
        if (opType == Filling || opType == Ident)
            m_spin = None;
        else
            m_spin = SymOp.getSpin();
        m_site = SymOp.getSite();
        m_opType = opType;
    }

    // Copy Constructor
    SymbolicOperator(const SymbolicOperator &SymOp) {
        m_part_type = SymOp.getPartType();
        m_spin = SymOp.getSpin();
        m_site = SymOp.getSite();
        m_opType = SymOp.getOpType();
    }
    // +---------+
    // | Getter |
    // +---------+
    m_SPIN getSpin() const { return m_spin; }
    pos_t getSite() const { return m_site; }
    m_OPTYPES getOpType() const { return m_opType; }
    part_type getPartType() const { return m_part_type; }

    // +---------+
    // | Methods |
    // +---------+
    /**! \brief
     *  This method prints the operator in an algebraic notation.
     */
    void print() const {
        std::string op_str;
        if (m_opType == Filling)
            op_str.append("F");
        else if (m_opType == Ident)
            op_str.append("I");
        else if (m_spin == Up || m_spin == Down)
            op_str.append("a");
        else if (m_spin == Zero)
            op_str.append("b");
        if (m_opType == Create)
            op_str.append("+");
        op_str.append("^");
        std::stringstream ss;
        ss << m_part_type;
        op_str.append(ss.str());
        op_str.append("_");
        ss.str(std::string());
        ss << m_site;
        op_str.append(ss.str());
        if (m_spin == Up)
            op_str.append("-up");
        else if (m_spin == Down)
            op_str.append("-down");
        std::cout << op_str;
    }
};


/**
 * @class SymbolicJordanWigner nu1_SymbolicJordanWigner.hpp
 * @brief This class handels the Jordan-Wigner transformation, based on the `SymbolicOperator`.
 * @note This is the most important and delicate part of the pre-BO model.
 * For a description of the algorithm, see my Master's thesis:
 * Correlation effects in Multicomponent Quantum Chemistry (2020), Robin Feldmann, Markus Reiher group ETH Zurich, p 27
 */
class SymbolicJordanWigner {

    typedef typename Lattice::part_type
            part_type; // unsigned integer denoting the particle type
    typedef typename Lattice::pos_t pos_t;

private:
    std::vector<SymbolicOperator> m_SymOpStr;
    std::vector<std::pair<part_type, pos_t>> m_nb_term;
    bool verbose = false;

public:
    // Default constructor
    SymbolicJordanWigner() {

    }
    // Constructor
    explicit SymbolicJordanWigner(std::vector<SymbolicOperator> SymOpStr) {

        if (SymOpStr.size()<1) throw std::runtime_error("Empty operator string in JW.");
        std::vector<part_type> vec_pt_check;
        // Identity check:
        // If the operator string contains only identities, now JW has to be applied.
        bool allIdent=true;
        for (auto const& Op : SymOpStr) {
            vec_pt_check.push_back(Op.getPartType());
            if (allIdent && Op.getOpType()!=SymbolicOperator::Ident) allIdent=false;
        }
        if (!(std::equal(vec_pt_check.begin() + 1, vec_pt_check.end(), vec_pt_check.begin())))
            throw std::runtime_error("JW must only be applied to operator strings containing equal particle types.");

        // Check passed.

        m_SymOpStr = SymOpStr;

        if (!allIdent) {

            // Identity removal
            for (auto it=m_SymOpStr.begin(); it!=m_SymOpStr.end();) {
                if (it->getOpType() == SymbolicOperator::Ident)
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

    /**! \brief
     * This getter returns the symbolic operator string
     * @return symbolic operator string
     */
    const std::vector<SymbolicOperator> &getSymOpStr() const { return m_SymOpStr; }

private:
    // +------------------------+
    // | Jordan--Wigner Filling |
    // +----------------------- +
    /*! \brief
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
            if (symOp.getSpin() == SymbolicOperator::Down ||
                symOp.getSpin() == SymbolicOperator::Up) {

                // c-(d) := F a(d)
                if (symOp.getSpin() == SymbolicOperator::Down &&
                    symOp.getOpType() == SymbolicOperator::Annihilate)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::Filling));

                if (lindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::Filling));

                new_SymOpStr.push_back(SymbolicOperator(symOp));

                // c+(d) := a+(d) F
                if (symOp.getSpin() == SymbolicOperator::Down &&
                    symOp.getOpType() == SymbolicOperator::Create)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::Filling));

                if (rindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::Filling));
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
                i->getOpType() == SymbolicOperator::Filling) {
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
