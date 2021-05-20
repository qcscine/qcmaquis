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


#include "nu1_SymbolicJordanWigner.hpp"

#include <map>


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
    std::vector<int> m_vec_orbitals; // Number of orbitals per type
    std::vector<std::pair<part_type, pos_t>> m_nbody_term;
    std::vector<size_t> m_spin; // physical spin is spin/2. So electrons have spin=1
    // If the input is in the Fiedler ordering:
    // std::vector<pos_t> m_order;
    // If the input is the standard ordering:
    std::vector<pos_t> m_inv_order;
    std::vector<std::vector<SymbolicOperator>> m_vec_SymOpStr;
    size_t m_NpartTypes;
    std::map<part_type, size_t> m_opPairMap;
    std::vector<std::vector<Spin>>
            m_spin_configs; // {SymbolicOperator::Down, SymbolicOperator::Up};
    std::vector<std::vector<std::vector<SymbolicOperator>>> m_allSymOpStr;

public:
    // +-------------+
    // | Constructor |
    // +-------------+
    NBodyTerm(std::vector<std::pair<part_type, pos_t>> nbody_term_,
              const std::vector<bool> &isFermion_, const std::vector<int>& vec_orbitals_,
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
        // std::pairs for each type. This requires a normal ordered input!!! N[ a+ b b+
        // a] = a+ b+ b a
        for (size_t i = 0; i < m_term_length / 2; i++)
            m_opPairMap[m_nbody_term[i].first] += 1;

        // Generate std::vector of nbody_terms only containing particles of one type
        std::vector<std::vector<std::pair<part_type, pos_t>>> vec_nbody_terms;
        for (auto const &iter : m_opPairMap) {
            std::vector<std::pair<part_type, pos_t>> temp_nbody;
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
            std::vector<std::vector<SymbolicOperator>> res_SymOpSts =
                    nBodyTerm2OperatorString(iter);
            m_allSymOpStr.push_back(res_SymOpSts);
        }

        // Create all combinations of spins between all particle types
        m_NpartTypes = vec_nbody_terms.size();
        int n = m_NpartTypes;
        std::vector<size_t> arr(n);
        generatePartConfigs(n, arr, 0);

    } // constructor

    // +--------+
    // | Getter |
    // +--------+
    const std::vector<std::vector<SymbolicOperator>> &getVecSymOpStr() const {
        return m_vec_SymOpStr;
    }

private:
    /*! \brief
     * This method takes the nbody term of a given particle type and creates the
     * strings of symbolic operators for all spin configurations. Then it also
     * takes care of the JW if the particle type is fermionic. It assumes that the
     * term is normal ordered.
     */
    std::vector<std::vector<SymbolicOperator>>
    nBodyTerm2OperatorString(std::vector<std::pair<part_type, pos_t>> temp_nbody) {
        size_t length = temp_nbody.size();
        part_type nt = temp_nbody[0].first;
        pos_t orb_index = 0;
        pos_t abs_index = 0;
        // -- Reset spin configs --
        m_spin_configs.resize(0);
        std::vector<std::vector<SymbolicOperator>> res_SymOpStr;
        // Boson S=0
        if (m_spin[nt] == 0) {
            std::vector<SymbolicOperator> tmp;
            for (size_t i = 0; i < temp_nbody.size(); i++) {
                orb_index = temp_nbody[i].second;
                abs_index = retrieve_abs_index(orb_index, nt);
                if (i < length / 2) {
                    SymbolicOperator tempOp(abs_index, OpType::Create, nt,
                                            Spin::Zero);
                    tmp.push_back(tempOp);
                }
                if (i >= length / 2) {
                    SymbolicOperator tempOp(abs_index, OpType::Annihilate, nt,
                                            Spin::Zero);
                    tmp.push_back(tempOp);
                }
            }
            res_SymOpStr.push_back(tmp);
        }
            // Fermion
        else if (m_spin[nt] == 1) {

            int n = length / 2;
            std::vector<size_t> arr(n);
            // This generates all possible Spin Configurations
            generateAllSpinConfigs(n, arr, 0);

            for (auto &spin_config : m_spin_configs) {
                for (int j = spin_config.size(); j-- != 0;) {
                    spin_config.push_back(spin_config.at(j));
                }
            }
            for (auto &spin_config : m_spin_configs) {
                std::vector<SymbolicOperator> tmp;
                for (size_t i = 0; i < temp_nbody.size(); i++) {
                    orb_index = temp_nbody[i].second;
                    abs_index = retrieve_abs_index(orb_index, nt);
                    // create
                    if (i < length / 2) {
                        SymbolicOperator tempOp(abs_index, OpType::Create, nt,
                                                spin_config.at(i));
                        tmp.push_back(tempOp);
                    }
                    if (i >= length / 2) {
                        SymbolicOperator tempOp(abs_index, OpType::Annihilate, nt,
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
     * std::vector of spin configurations.
     */
    void AddToConfigs(std::vector<size_t> arr, int n) {
        std::vector<Spin> temp_config;
        for (int i = 0; i < n; i++) {
            if (arr[i] == 0)
                temp_config.push_back(Spin::Down);
            else if (arr[i] == 1)
                temp_config.push_back(Spin::Up);
        }
        m_spin_configs.push_back(temp_config);
    }

    /*! \brief
     * This method recursively constructs all possible spin configurations for a
     * given particle type.
     */
    void generateAllSpinConfigs(int n, std::vector<size_t> arr, int i) {
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
    void AddToPartConfigs(std::vector<size_t> arr, int n) {
        std::vector<SymbolicOperator> temp_config;

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
    void generatePartConfigs(int n, std::vector<size_t> arr, int i) {
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
    pos_t retrieve_abs_index(const pos_t& rel_orb_index, const part_type& type) {
        pos_t abs_index = 0;
        for (size_t i = 0; i < type; i++) {
            abs_index += m_vec_orbitals[i];
        }
        abs_index += rel_orb_index;
        return m_inv_order[abs_index];
    }

public:
    static pos_t retrieve_abs_index(const pos_t& rel_orb_index, const part_type& type, const std::vector<int>& vec_orbitals,
                                    const std::vector<pos_t>& inv_order) {
        pos_t abs_index = 0;
        for (size_t i = 0; i < type; i++) {
            abs_index += vec_orbitals[i];
        }
        abs_index += rel_orb_index;
        return inv_order[abs_index];
    }


};

#endif // MAQUIS_DMRG_NU1_NBODYTERM_HPP
