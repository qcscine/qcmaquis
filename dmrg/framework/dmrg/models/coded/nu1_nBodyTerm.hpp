/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2017 by Alberto Baiardi <alberto.baiardi@phys.chem.ethz.ch>
 *               2020 by Robin Feldmann <robinfe@student.ethz.ch>
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

#ifndef MAQUIS_DMRG_NU1_NBODYTERM_HPP
#define MAQUIS_DMRG_NU1_NBODYTERM_HPP

#include <dmrg/models/lattice.h>

#include <map>

using std::pair;
using std::vector;

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
    enum m_SPIN { UP, DOWN, ZERO, NONE };
    enum m_OPTYPES { CREATE, ANNIHILATE, FILLING, IDENT };

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
                     m_SPIN spin = NONE)
            : m_site(site), m_opType(opType), m_part_type(pt), m_spin(spin) {}

    // Constructor
    SymbolicOperator(const SymbolicOperator &SymOp, m_OPTYPES opType) {
        m_part_type = SymOp.getPartType();
        if (opType == FILLING || opType == IDENT)
            m_spin = NONE;
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
        if (m_opType == FILLING)
            op_str.append("F");
        else if (m_opType == IDENT)
            op_str.append("I");
        else if (m_spin == UP || m_spin == DOWN)
            op_str.append("a");
        else if (m_spin == ZERO)
            op_str.append("b");
        if (m_opType == CREATE)
            op_str.append("+");
        op_str.append("^");
        std::stringstream ss;
        ss << m_part_type;
        op_str.append(ss.str());
        op_str.append("_");
        ss.str(std::string());
        ss << m_site;
        op_str.append(ss.str());
        if (m_spin == UP)
            op_str.append("-up");
        else if (m_spin == DOWN)
            op_str.append("-down");
        std::cout << op_str;
    }
};

/**
 * @class SymbolicJordanWigner nu1_nBodyterm.hpp
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
    vector<SymbolicOperator> m_SymOpStr;
    vector<pair<part_type, pos_t>> m_nb_term;
    bool verbose = false;

public:
    // Default constructor
    SymbolicJordanWigner() {

    }
    // Constructor
    explicit SymbolicJordanWigner(vector<SymbolicOperator> SymOpStr) {

        if (SymOpStr.size()<1) throw std::runtime_error("Empty operator string in JW.");
        vector<part_type> vec_pt_check;
        // Identity check:
        // If the operator string contains only identities, now JW has to be applied.
        bool allIdent=true;
        for (auto const& Op : SymOpStr) {
            vec_pt_check.push_back(Op.getPartType());
            if (allIdent && Op.getOpType()!=SymbolicOperator::IDENT) allIdent=false;
        }
        if (!(std::equal(vec_pt_check.begin() + 1, vec_pt_check.end(), vec_pt_check.begin())))
            throw std::runtime_error("JW must only be applied to operator strings containing equal particle types.");

        // Check passed.

        m_SymOpStr = SymOpStr;

        if (!allIdent) {

            // Identity removal
            for (auto it=m_SymOpStr.begin(); it!=m_SymOpStr.end();) {
                if (it->getOpType() == SymbolicOperator::IDENT)
                    it = m_SymOpStr.erase(it);
                else
                    ++it;
            }
            for (auto const& iter : m_SymOpStr)
                m_nb_term.push_back(std::make_pair(iter.getPartType(), iter.getSite()));

            vector<pos_t> lindex;
            vector<pos_t> rindex;

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
    const vector<SymbolicOperator> &getSymOpStr() const { return m_SymOpStr; }

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
    JW_generate_lr_index(const vector<pair<part_type, pos_t>> &nb_term,
                         vector<pos_t> &lindex, vector<pos_t> &rindex) {
        unsigned int length = nb_term.size();
        for (unsigned int i = 0; i < length; i++) {
            lindex.push_back(0);
            rindex.push_back(0);
        }
        //
        // VERY IMPORTANT:  This vector of maps takes care of the Jordan--Wigner
        // transformation!!
        //                  --> the procedure is very delicate and error-prone,
        //                  --> any modifications must be debugged very carefully
        //
        // This vector of maps stores the filling string for the relevant sites for
        // each operator. key: index, value: true if operator at position in vector
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

    void insertFillings(const vector<pos_t> &lindex,
                        const vector<pos_t> &rindex) {
        vector<SymbolicOperator> new_SymOpStr;
        size_t count = 0;
        for (auto const &symOp : m_SymOpStr) {
            // fermion
            if (symOp.getSpin() == SymbolicOperator::DOWN ||
                symOp.getSpin() == SymbolicOperator::UP) {

                // c-(d) := F a(d)
                if (symOp.getSpin() == SymbolicOperator::DOWN &&
                    symOp.getOpType() == SymbolicOperator::ANNIHILATE)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::FILLING));

                if (lindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::FILLING));

                new_SymOpStr.push_back(SymbolicOperator(symOp));

                // c+(d) := a+(d) F
                if (symOp.getSpin() == SymbolicOperator::DOWN &&
                    symOp.getOpType() == SymbolicOperator::CREATE)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::FILLING));

                if (rindex[count] == 1)
                    new_SymOpStr.push_back(
                            SymbolicOperator(symOp, SymbolicOperator::FILLING));
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
                i->getOpType() == SymbolicOperator::FILLING) {
                i = new_SymOpStr.erase(i);
                i = new_SymOpStr.erase(i);
            } else {
                i++;
            }
        }
        m_SymOpStr = new_SymOpStr;
    }
};

/**
 * @class NBodyTerm nu1_nBodyTerm.hpp
 * @brief Takes an n-body 2nd quantization string that can contain operators of different particle types and creates
 *        a string of `SymbolicOperator`s. Moreover, it constructs all different spin-combinations.
 *        Note that creation operators are always on the left, and -obviously- there must be even numbers of operators
 *        for all types.
 */
class NBodyTerm {
    // Types definition
    typedef typename Lattice::part_type part_type; // unsigned integer denoting the particle type
    typedef typename Lattice::pos_t pos_t;

private:
    // +---------+
    // | Members |
    // +---------+
    vector<unsigned int> m_vec_orbitals; // Number of orbitals per type
    vector<pair<part_type, pos_t>> m_nbody_term;
    vector<size_t> m_spin; // physical spin is spin/2. So electrons have spin=1
    // If the input is in the Fiedler ordering:
    // std::vector<pos_t> m_order;
    // If the input is the standard ordering:
    std::vector<pos_t> m_inv_order;
    vector<vector<SymbolicOperator>> m_vec_SymOpStr;
    size_t m_NpartTypes;
    std::map<part_type, size_t> m_opPairMap;
    vector<vector<SymbolicOperator::m_SPIN>>
            m_spin_configs; // {SymbolicOperator::DOWN, SymbolicOperator::UP};
    vector<vector<vector<SymbolicOperator>>> m_allSymOpStr;

public:
    // +-------------+
    // | Constructor |
    // +-------------+
    NBodyTerm(vector<pair<part_type, pos_t>> nbody_term_,
              const vector<bool> &isFermion_, vector<unsigned int> &vec_orbitals_,
              const std::vector<pos_t> &order_) {
        m_nbody_term = nbody_term_;
        size_t m_term_length = m_nbody_term.size();
        m_vec_orbitals = vec_orbitals_;
        m_inv_order = order_;

        // This is redundant right now.
        for (auto const &sigma : isFermion_)
            m_spin.push_back(sigma);

        // nbody term must have length%2=0!
        assert(m_term_length % 2 == 0);
        // Count the number of different particle types and the number of operator
        // pairs for each type. This requires a normal ordered input!!! N[ a+ b b+
        // a] = a+ b+ b a
        for (size_t i = 0; i < m_term_length / 2; i++)
            m_opPairMap[m_nbody_term[i].first] += 1;

        // Generate vector of nbody_terms only containing particles of one type
        vector<vector<pair<part_type, pos_t>>> vec_nbody_terms;
        for (auto const &iter : m_opPairMap) {
            vector<pair<part_type, pos_t>> temp_nbody;
            for (auto i = m_nbody_term.begin(); i != m_nbody_term.end();) {
                // auto n = std::next(i);
                if (i == m_nbody_term.end())
                    break;
                if (iter.first == i->first) {
                    temp_nbody.push_back(*i);
                    i = m_nbody_term.erase(i);
                } else {
                    i++;
                }
            }
            vec_nbody_terms.push_back(temp_nbody);
        }

        // Generate all operator strings with spin configurations for all particle
        // types
        for (auto const &iter : vec_nbody_terms) {
            vector<vector<SymbolicOperator>> res_SymOpSts =
                    nBodyTerm2OperatorString(iter);
            m_allSymOpStr.push_back(res_SymOpSts);
        }

        // Create all combinations of spins between all particle types
        m_NpartTypes = vec_nbody_terms.size();
        int n = m_NpartTypes;
        vector<size_t> arr(n);
        generatePartConfigs(n, arr, 0);

    } // constructor

    // +--------+
    // | Getter |
    // +--------+
    const vector<vector<SymbolicOperator>> &getVecSymOpStr() const {
        return m_vec_SymOpStr;
    }

private:
    /*! \brief
     * This method takes the nbody term of a given particle type and creates the
     * strings of symbolic operators for all spin configurations. Then it also
     * takes care of the JW if the particle type is fermionic. It assumes that the
     * term is normal ordered.
     */
    vector<vector<SymbolicOperator>>
    nBodyTerm2OperatorString(vector<pair<part_type, pos_t>> temp_nbody) {
        size_t length = temp_nbody.size();
        part_type nt = temp_nbody[0].first;
        pos_t orb_index = 0;
        pos_t abs_index = 0;
        // -- Reset spin configs --
        m_spin_configs.resize(0);
        vector<vector<SymbolicOperator>> res_SymOpStr;
        // Boson S=0
        if (m_spin[nt] == 0) {
            vector<SymbolicOperator> tmp;
            for (size_t i = 0; i < temp_nbody.size(); i++) {
                orb_index = temp_nbody[i].second;
                abs_index = retrieve_abs_index(orb_index, nt);
                if (i < length / 2) {
                    SymbolicOperator tempOp(abs_index, SymbolicOperator::CREATE, nt,
                                            SymbolicOperator::ZERO);
                    tmp.push_back(tempOp);
                }
                if (i >= length / 2) {
                    SymbolicOperator tempOp(abs_index, SymbolicOperator::ANNIHILATE, nt,
                                            SymbolicOperator::ZERO);
                    tmp.push_back(tempOp);
                }
            }
            res_SymOpStr.push_back(tmp);
        }
            // Fermion
        else if (m_spin[nt] == 1) {

            int n = length / 2;
            vector<size_t> arr(n);
            // This generates all possible Spin Configurations
            generateAllSpinConfigs(n, arr, 0);

            for (auto &spin_config : m_spin_configs) {
                for (unsigned j = spin_config.size(); j-- != 0;) {
                    spin_config.push_back(spin_config.at(j));
                }
            }
            for (auto &spin_config : m_spin_configs) {
                vector<SymbolicOperator> tmp;
                for (size_t i = 0; i < temp_nbody.size(); i++) {
                    orb_index = temp_nbody[i].second;
                    abs_index = retrieve_abs_index(orb_index, nt);
                    // create
                    if (i < length / 2) {
                        SymbolicOperator tempOp(abs_index, SymbolicOperator::CREATE, nt,
                                                spin_config.at(i));
                        tmp.push_back(tempOp);
                    }
                    if (i >= length / 2) {
                        SymbolicOperator tempOp(abs_index, SymbolicOperator::ANNIHILATE, nt,
                                                spin_config.at(i));
                        tmp.push_back(tempOp);
                    }
                }
                SymbolicJordanWigner JW(tmp);
                res_SymOpStr.push_back(JW.getSymOpStr());
            }
        }
        return res_SymOpStr;
    }

    /*! \brief
     * This method takes the indices of the permutation and builds the according
     * vector of spin configurations.
     */
    void AddToConfigs(vector<size_t> arr, int n) {
        vector<SymbolicOperator::m_SPIN> temp_config;
        for (int i = 0; i < n; i++) {
            if (arr[i] == 0)
                temp_config.push_back(SymbolicOperator::DOWN);
            else if (arr[i] == 1)
                temp_config.push_back(SymbolicOperator::UP);
        }
        m_spin_configs.push_back(temp_config);
    }

    /*! \brief
     * This method recursively constructs all possible spin configurations for a
     * given particle type.
     */
    void generateAllSpinConfigs(int n, vector<size_t> arr, int i) {
        if (i == n) {
            AddToConfigs(arr, n);
            return;
        }
        // First assign "0" at ith position
        // and try for all other permutations
        // for remaining positions
        arr[i] = 0;
        generateAllSpinConfigs(n, arr, i + 1);
        // And then assign "1" at ith position
        // and try for all other permutations
        // for remaining positions
        arr[i] = 1;
        generateAllSpinConfigs(n, arr, i + 1);
    }

    /*! \brief
     * This method takes the indices of the permutation and builds the according
     * string of symbolic operators.
     */
    void AddToPartConfigs(vector<size_t> arr, int n) {
        vector<SymbolicOperator> temp_config;

        for (int i = 0; i < n; i++) {
            for (const SymbolicOperator &SymOp : m_allSymOpStr.at(i).at(arr[i]))
                temp_config.push_back(SymOp);
        }

        m_vec_SymOpStr.push_back(temp_config);
    }
    /*! \brief
     * This method recursively constructs all possible configurations how the spin
     * configurations of different particle types can be combined.
     */
    void generatePartConfigs(int n, vector<size_t> arr, int i) {
        if (i == n) {
            AddToPartConfigs(arr, n);
            return;
        }
        for (size_t count = 0; count < m_allSymOpStr.at(i).size(); count++) {
            arr[i] = count;
            generatePartConfigs(n, arr, i + 1);
        }
    }

    // +-------------------------+
    // | Retrieve Absolute Index |
    // +-------------------------+
    /*! \brief
     * This method uses the type and the orbital index to retrieve the index of
     * the orbital on the lattice. If the lattice is permuted the function takes
     * care of it.
     */
    pos_t retrieve_abs_index(pos_t rel_orb_index, part_type type) {
        pos_t abs_index = 0;
        for (size_t i = 0; i < type; i++) {
            abs_index += m_vec_orbitals[i];
        }
        abs_index += rel_orb_index;
        return m_inv_order[abs_index];
    }
};

#endif // MAQUIS_DMRG_NU1_NBODYTERM_HPP
