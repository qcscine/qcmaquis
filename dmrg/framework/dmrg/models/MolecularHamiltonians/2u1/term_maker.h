/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef QC_TERMMAKER_H
#define QC_TERMMAKER_H

#include "dmrg/models/OperatorHandlers/TagHandler.h"
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/models/term_descriptor.h"

template <class M, class S>
struct TermMaker {
  // == Types declaration ==
  using pos_t = typename Lattice::pos_t;
  using value_type = typename M::value_type;
  using term_descriptor = ::term_descriptor<value_type>;
  using tag_type = typename TagHandler<M, S>::tag_type;
  using pos_op_t = typename term_descriptor::value_type;
  using sc_t = typename S::subcharge;

  /**
   * @brief Generates the non-equivalent indexes for a given one-electron
   * integral.
   */
  static std::vector<std::array<int, 2> > generateTwofoldSymmetricIndex(
      int i, int j
  ) {
    using RetType = std::array<int, 2>;
    RetType ret = {i, j};
    // Generates the permutation
    std::array<int, 2> ret1 = {ret[0], ret[1]}, ret2 = {ret[1], ret[0]};
    std::set<std::array<int, 2> > tmp = {ret1, ret2};
    return std::vector<RetType>(tmp.begin(), tmp.end());
  }

  /**
   * @brief Generates the non-equivalent indexes for a given two-electron
   * integral. Note that the input is assumed to be given in chemical notation,
   * i.e. (ij,kl) == <ik||jl>.
   */
  static std::vector<std::array<int, 4> > generateEightfoldSymmetricIndex(
      int i, int j, int k, int l
  ) {
    using RetType = std::array<int, 4>;
    RetType ret = {i, j, k, l};
    // Generates the permutation
    std::array<int, 4> ret1 = {ret[0], ret[1], ret[3], ret[2]},
                       ret2 = {ret[1], ret[0], ret[2], ret[3]},
                       ret3 = {ret[1], ret[0], ret[3], ret[2]},
                       ret4 = {ret[2], ret[3], ret[0], ret[1]},
                       ret5 = {ret[2], ret[3], ret[1], ret[0]},
                       ret6 = {ret[3], ret[2], ret[0], ret[1]},
                       ret7 = {ret[3], ret[2], ret[1], ret[0]};
    std::set<std::array<int, 4> > tmp = {ret,  ret1, ret2, ret3,
                                         ret4, ret5, ret6, ret7};
    return std::vector<RetType>(tmp.begin(), tmp.end());
  }

  /**
   * @brief Generates the non-equivalent indexes for a given three-electron
   * integral. Note that the input is assumed to be given in chemical notation,
   * i.e. (ij,kl,mn), where ij --> r1, kl --> r2, mn --> r3
   */
  static std::vector<std::array<int, 6> > generateThreeBodySymmetricIndex(
      int i, int j, int k, int l, int m, int n
  ) {
    using RetType = std::array<int, 6>;
    RetType ret = {i, j, k, l, m, n};
    // Generates the permutation
    std::vector<std::array<int, 6> > terms = {
        // 123
        {ret[0], ret[1], ret[2], ret[3], ret[4], ret[5]},
        {ret[0], ret[1], ret[2], ret[3], ret[5], ret[4]},
        {ret[0], ret[1], ret[3], ret[2], ret[4], ret[5]},
        {ret[0], ret[1], ret[3], ret[2], ret[5], ret[4]},
        {ret[1], ret[0], ret[2], ret[3], ret[4], ret[5]},
        {ret[1], ret[0], ret[2], ret[3], ret[5], ret[4]},
        {ret[1], ret[0], ret[3], ret[2], ret[4], ret[5]},
        {ret[1], ret[0], ret[3], ret[2], ret[5], ret[4]},
        // 213
        {ret[2], ret[3], ret[0], ret[1], ret[4], ret[5]},
        {ret[2], ret[3], ret[0], ret[1], ret[5], ret[4]},
        {ret[3], ret[2], ret[0], ret[1], ret[4], ret[5]},
        {ret[3], ret[2], ret[0], ret[1], ret[5], ret[4]},
        {ret[2], ret[3], ret[1], ret[0], ret[4], ret[5]},
        {ret[2], ret[3], ret[1], ret[0], ret[5], ret[4]},
        {ret[3], ret[2], ret[1], ret[0], ret[4], ret[5]},
        {ret[3], ret[2], ret[1], ret[0], ret[5], ret[4]},
        // 132
        {ret[0], ret[1], ret[4], ret[5], ret[2], ret[3]},
        {ret[0], ret[1], ret[5], ret[4], ret[2], ret[3]},
        {ret[0], ret[1], ret[4], ret[5], ret[3], ret[2]},
        {ret[0], ret[1], ret[5], ret[4], ret[3], ret[2]},
        {ret[1], ret[0], ret[4], ret[5], ret[2], ret[3]},
        {ret[1], ret[0], ret[5], ret[4], ret[2], ret[3]},
        {ret[1], ret[0], ret[4], ret[5], ret[3], ret[2]},
        {ret[1], ret[0], ret[5], ret[4], ret[3], ret[2]},
        // 231
        {ret[2], ret[3], ret[4], ret[5], ret[0], ret[1]},
        {ret[2], ret[3], ret[5], ret[4], ret[0], ret[1]},
        {ret[3], ret[2], ret[4], ret[5], ret[0], ret[1]},
        {ret[3], ret[2], ret[5], ret[4], ret[0], ret[1]},
        {ret[2], ret[3], ret[4], ret[5], ret[1], ret[0]},
        {ret[2], ret[3], ret[5], ret[4], ret[1], ret[0]},
        {ret[3], ret[2], ret[4], ret[5], ret[1], ret[0]},
        {ret[3], ret[2], ret[5], ret[4], ret[1], ret[0]},
        // 312
        {ret[4], ret[5], ret[0], ret[1], ret[2], ret[3]},
        {ret[5], ret[4], ret[0], ret[1], ret[2], ret[3]},
        {ret[4], ret[5], ret[0], ret[1], ret[3], ret[2]},
        {ret[5], ret[4], ret[0], ret[1], ret[3], ret[2]},
        {ret[4], ret[5], ret[1], ret[0], ret[2], ret[3]},
        {ret[5], ret[4], ret[1], ret[0], ret[2], ret[3]},
        {ret[4], ret[5], ret[1], ret[0], ret[3], ret[2]},
        {ret[5], ret[4], ret[1], ret[0], ret[3], ret[2]},
        // 321
        {ret[4], ret[5], ret[2], ret[3], ret[0], ret[1]},
        {ret[5], ret[4], ret[2], ret[3], ret[0], ret[1]},
        {ret[4], ret[5], ret[3], ret[2], ret[0], ret[1]},
        {ret[5], ret[4], ret[3], ret[2], ret[0], ret[1]},
        {ret[4], ret[5], ret[2], ret[3], ret[1], ret[0]},
        {ret[5], ret[4], ret[2], ret[3], ret[1], ret[0]},
        {ret[4], ret[5], ret[3], ret[2], ret[1], ret[0]},
        {ret[5], ret[4], ret[3], ret[2], ret[1], ret[0]}};
    std::set<std::array<int, 6> > tmp;
    for (const auto& iTerm : terms) tmp.insert(iTerm);
    return std::vector<RetType>(tmp.begin(), tmp.end());
  }

  /**
   * @brief Generates the non-equivalent indexes based on two-fold symmetry.
   * The only symmetry exploited is (ij|kl) = (kl|ij).
   */
  static std::vector<std::array<int, 4> > generateTwofoldSymmetricIndex(
      int i, int j, int k, int l
  ) {
    using RetType = std::array<int, 4>;
    RetType ret = {i, j, k, l};
    // Generates the permutation
    std::array<int, 4> ret1 = {ret[2], ret[3], ret[0], ret[1]};
    std::set<std::array<int, 4> > tmp = {ret, ret1};
    return std::vector<RetType>(tmp.begin(), tmp.end());
  }

  /** @brief Comparison operator between two tags */
  static bool compare_tag(pos_op_t p1, pos_op_t p2) {
    return std::get<0>(p1) < std::get<0>(p2);
  }

  /** @brief Function to add a simple two-term operator */
  static term_descriptor two_term(
      bool sign, std::vector<tag_type> const& fill_op, value_type scale,
      pos_t i, pos_t j, std::vector<tag_type> const& op1,
      std::vector<tag_type> const& op2,
      std::shared_ptr<TagHandler<M, S> > op_table, Lattice const& lat
  ) {
    term_descriptor term;
    term.is_fermionic = sign;
    term.coeff = scale;
    term.push_back(std::make_pair(i, op1[lat.get_prop<sc_t>("type", i)]));
    term.push_back(std::make_pair(j, op2[lat.get_prop<sc_t>("type", j)]));
    return term;
  }

  /**
   * @brief Same as above, but for positional operators.
   * Note that here, by positional, we mean that also the filling operator is
   * added in the operator evaluation, to take the Jordan-Wigner transformation
   * into account.
   */
  static term_descriptor positional_two_term(
      bool sign, std::vector<tag_type> const& fill_op, value_type scale,
      pos_t i, pos_t j, std::vector<tag_type> const& op1,
      std::vector<tag_type> const& op2,
      std::shared_ptr<TagHandler<M, S> > op_table, Lattice const& lat
  ) {
    term_descriptor term;
    term.is_fermionic = sign;
    term.coeff = scale;

    std::pair<tag_type, value_type> ptag;
    if (i < j) {
      ptag = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", i)],
          op1[lat.get_prop<sc_t>("type", i)]
      );
      term.push_back(std::make_pair(i, ptag.first));
      term.push_back(std::make_pair(j, op2[lat.get_prop<sc_t>("type", j)]));
      term.coeff *= ptag.second;
    } else {
      ptag = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", j)],
          op2[lat.get_prop<sc_t>("type", j)]
      );
      term.push_back(std::make_pair(i, op1[lat.get_prop<sc_t>("type", i)]));
      term.push_back(std::make_pair(j, ptag.first));
      term.coeff *= -ptag.second;
    }
    return term;
  }

  /** @brief Same as above, but has three operators, where the first two act on
   * the same site and can be multiplied together */
  static term_descriptor positional_two_term(
      bool sign, std::vector<tag_type> const& fill_op, value_type scale,
      pos_t i, pos_t j, std::vector<tag_type> const& op1,
      const std::vector<tag_type>& op2, const std::vector<tag_type>& op3,
      std::shared_ptr<TagHandler<M, S> > op_table, const Lattice& lat
  ) {
    term_descriptor term;
    term.is_fermionic = sign;
    term.coeff = scale;
    std::pair<tag_type, value_type> pre_ptag;
    pre_ptag = op_table->get_product_tag(
        op1[lat.get_prop<sc_t>("type", i)], op2[lat.get_prop<sc_t>("type", i)]
    );
    std::pair<tag_type, value_type> ptag;
    if (i < j) {
      ptag = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", i)], pre_ptag.first
      );
      term.push_back(std::make_pair(i, ptag.first));
      term.push_back(std::make_pair(j, op3[lat.get_prop<sc_t>("type", j)]));
      term.coeff *= ptag.second * pre_ptag.second;
    } else {
      ptag = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", j)],
          op3[lat.get_prop<sc_t>("type", j)]
      );
      term.push_back(std::make_pair(i, pre_ptag.first));
      term.push_back(std::make_pair(j, ptag.first));
      term.coeff *= -ptag.second * pre_ptag.second;
    }
    return term;
  }

  /** @brief Three-site operator. Here the filling is added a-priori */
  static term_descriptor three_term(
      std::vector<tag_type> const& ident, std::vector<tag_type> const& fill_op,
      value_type scale, pos_t pb, pos_t p1, pos_t p2,
      std::vector<tag_type> const& opb1, std::vector<tag_type> const& opb2,
      std::vector<tag_type> const& ops1, std::vector<tag_type> const& ops2,
      std::shared_ptr<TagHandler<M, S> > op_table, Lattice const& lat
  ) {
    // Set up the descriptor object
    term_descriptor term;
    term.is_fermionic = true;
    term.coeff = scale;

    tag_type boson_op;
    tag_type op1 = ops1[lat.get_prop<sc_t>("type", p1)];
    tag_type op2 = ops2[lat.get_prop<sc_t>("type", p2)];
    std::pair<tag_type, value_type> ptag1, ptag2;

    if ((pb > p1 && pb < p2) || (pb > p2 && pb < p1)) {
      // if the bosonic operator is in between
      // the fermionic operators, multiply with fill
      ptag1 = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", pb)],
          opb2[lat.get_prop<sc_t>("type", pb)]
      );
      term.coeff *= ptag1.second;
      ptag2 = op_table->get_product_tag(
          ptag1.first, opb1[lat.get_prop<sc_t>("type", pb)]
      );
      term.coeff *= ptag2.second;
      boson_op = ptag2.first;
    } else {
      ptag1 = op_table->get_product_tag(
          opb2[lat.get_prop<sc_t>("type", pb)],
          opb1[lat.get_prop<sc_t>("type", pb)]
      );
      boson_op = ptag1.first;
      term.coeff *= ptag1.second;
    }

    if (p1 < p2) {
      ptag1 = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", p1)],
          ops1[lat.get_prop<sc_t>("type", p1)]
      );
      op1 = ptag1.first;
      term.coeff *= ptag1.second;
    } else {
      ptag1 = op_table->get_product_tag(
          fill_op[lat.get_prop<sc_t>("type", p2)],
          ops2[lat.get_prop<sc_t>("type", p2)]
      );
      op2 = ptag1.first;
      term.coeff *= -ptag1.second;
    }

    std::vector<pos_op_t> sterm;
    sterm.push_back(std::make_pair(pb, boson_op));
    sterm.push_back(std::make_pair(p1, op1));
    sterm.push_back(std::make_pair(p2, op2));
    std::sort(sterm.begin(), sterm.end(), compare_tag);

    term.push_back(sterm[0]);
    term.push_back(sterm[1]);
    term.push_back(sterm[2]);

    return term;
  }

  /** @brief Three-site operator. Here the filling is added a-priori */
  static term_descriptor four_term(
      std::vector<tag_type> const& ident, std::vector<tag_type> const& fill_op,
      value_type scale, pos_t i, pos_t j, pos_t k, pos_t l,
      std::vector<tag_type> const& op_i, std::vector<tag_type> const& op_j,
      std::vector<tag_type> const& op_k, std::vector<tag_type> const& op_l,
      std::shared_ptr<TagHandler<M, S> > op_table, Lattice const& lat
  ) {
    term_descriptor term;
    term.is_fermionic = true;
    term.coeff = scale;

    // Simple O(n^2) algorithm to determine sign of permutation
    pos_t idx[] = {i, j, k, l};
    pos_t inv_count = 0, n = 4;
    for (pos_t c1 = 0; c1 < n - 1; c1++)
      for (pos_t c2 = c1 + 1; c2 < n; c2++)
        if (idx[c1] > idx[c2]) inv_count++;

    std::vector<pos_op_t> sterm;
    sterm.push_back(std::make_pair(i, op_i[lat.get_prop<sc_t>("type", i)]));
    sterm.push_back(std::make_pair(j, op_j[lat.get_prop<sc_t>("type", j)]));
    sterm.push_back(std::make_pair(k, op_k[lat.get_prop<sc_t>("type", k)]));
    sterm.push_back(std::make_pair(l, op_l[lat.get_prop<sc_t>("type", l)]));
    std::sort(sterm.begin(), sterm.end(), compare_tag);

    std::pair<tag_type, value_type> ptag;
    ptag = op_table->get_product_tag(
        fill_op[lat.get_prop<sc_t>("type", std::get<0>(sterm[0]))],
        std::get<1>(sterm[0])
    );
    std::get<1>(sterm[0]) = ptag.first;
    term.coeff *= ptag.second;
    ptag = op_table->get_product_tag(
        fill_op[lat.get_prop<sc_t>("type", std::get<0>(sterm[2]))],
        std::get<1>(sterm[2])
    );
    std::get<1>(sterm[2]) = ptag.first;
    term.coeff *= ptag.second;

    if (inv_count % 2) term.coeff = -term.coeff;

    term.push_back(sterm[0]);
    term.push_back(sterm[1]);
    term.push_back(sterm[2]);
    term.push_back(sterm[3]);
    return term;
  }
};

#endif
