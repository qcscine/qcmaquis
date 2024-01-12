/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied
 * Biosciences, Reiher Group. See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_SQUARE_LATTICE
#define MAQUIS_DMRG_SQUARE_LATTICE

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"

class SquareLattice : public lattice_impl {
 public:
  SquareLattice(BaseParameters& parms)
      : L_(parms["L"]), W_(parms["W"]), a(parms["a"]) {}

  /*
   0 4  8 12
   1 5  9 13
   2 6 10 14
   3 7 11 15
   */
  std::vector<int> forward(int p) const override {
    std::vector<int> ret;
    if (p + 1 < L_ * W_ && (p + 1) % W_ != 0) {
      ret.push_back(p + 1);
    }
    if (p + W_ < L_ * W_) {
      ret.push_back(p + W_);
    }
    return ret;
  }

  std::vector<int> all(int p) const override {
    std::vector<int> ret = forward(p);
    if (p >= 1 && p % W_ != 0) {
      ret.push_back(p - 1);
    }
    if (p >= W_) {
      ret.push_back(p - W_);
    }
    return ret;
  }

  /** @brief Getter for the lattice size */
  int size() const override { return L_ * W_; }

  /** @brief Getter for the number of types of sites */
  int getMaxType() const override { return 1; }

  std::any get_prop_(std::string const& property, std::vector<pos_t> const& pos)
      const override {
    if (property == "label" && pos.size() == 1) {
      return std::any(site_label(pos[0]));
    } else if (property == "label" && pos.size() == 2) {
      return std::any(bond_label(pos[0], pos[1]));
    } else if (property == "type" && pos.size() == 1) {
      return std::any(0);
    } else if (property == "type" && pos.size() == 2) {
      return std::any(0);
    } else if (property == "x" && pos.size() == 1) {
      return std::any(x(pos[0]));
    } else if (property == "y" && pos.size() == 1) {
      return std::any(y(pos[0]));
    } else if (property == "wraps_pbc" && pos.size() == 2) {
      return std::any(false);
    } else if (property == "NumTypes") {
      return std::any(1);
    } else if (property == "ParticleType") {
      return std::any(0);
    } else {
      std::ostringstream ss;
      ss << "No property '" << property << "' with " << pos.size()
         << " points implemented.";
      throw std::runtime_error(ss.str());
      return std::any();
    }
  }

 private:
  double x(int i) const { return a * int(i / W_); }

  double y(int i) const { return a * (i % W_); }

  std::string site_label(int i) const {
    return "( " +
           (boost::lexical_cast<std::string>(x(i)) + "," +
            boost::lexical_cast<std::string>(y(i))) +
           " )";
  }

  std::string bond_label(int i, int j) const {
    return (
        "( " +
        (boost::lexical_cast<std::string>(x(i)) + "," +
         boost::lexical_cast<std::string>(y(i))) +
        " )" + " -- " + "( " +
        (boost::lexical_cast<std::string>(x(j)) + "," +
         boost::lexical_cast<std::string>(y(j))) +
        " )"
    );
  }

  int L_, W_;
  double a;
};

#endif
