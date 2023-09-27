/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_lattice_hpp
#define MAQUIS_DMRG_lattice_hpp

#include <numeric>
#include <sstream>
#include <set>
#include <vector>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include "dmrg/models/lattice/lattice.h"
#include "dmrg/utils/BaseParameters.h"

class ChainLattice : public lattice_impl
{
public:
  // Types definition
  using pos_t = typename lattice_impl::pos_t;

  /** @brief Constructor from a parameter container */
  ChainLattice (BaseParameters & parms, bool pbc_=false)
    : L(parms["L"]) , pbc(pbc_) { }

  /** @brief Constructor from parameters */
  ChainLattice (int L_, bool pbc_=false)
    : L(L_), pbc(pbc_) { }

  /** @brief Get next position in the lattice */
  std::vector<pos_t> forward(pos_t i) const
  {
    std::vector<pos_t> ret;
    if (i < L-1)
      ret.push_back(i+1);
    if (pbc && i == L-1)
      ret.push_back(0);
    return ret;
  }

  /** @brief Getter for the neighbours of a given site */
  std::vector<pos_t> all(pos_t i) const
  {
    std::vector<pos_t> ret;
    if (i < L-1)
      ret.push_back(i+1);
    if (i > 0)
      ret.push_back(i-1);
    if (pbc && i == L-1)
      ret.push_back(0);
    if (pbc && i == 0)
      ret.push_back(L-1);
    return ret;
  }

  /** @brief Getter for a generic property of the lattice */
  boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
  {
    if (property == "label" && pos.size() == 1)
      return boost::any( site_label(pos[0]) );
    else if (property == "label" && pos.size() == 2)
      return boost::any( bond_label(pos[0], pos[1]) );
    else if (property == "type" && pos.size() == 1)
      return boost::any( 0 );
    else if (property == "type" && pos.size() == 2)
      return boost::any( 0 );
    else if (property == "x" && pos.size() == 1)
      return boost::any( pos[0] );
    else if (property == "at_open_boundary" && pos.size() == 1)
      return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L-1) );
    else if (property == "at_open_left_boundary" && pos.size() == 1)
      return boost::any( (!pbc) && pos[0]==0 );
    else if (property == "at_open_right_boundary" && pos.size() == 1)
      return boost::any( (!pbc) && pos[0]==L-1 );
    else if (property == "wraps_pbc" && pos.size() == 2)
      return boost::any( (pos[0] < pos[1]) );
    else if (property == "NumTypes")
      return boost::any( 1 );
    else if (property == "ParticleType" && pos.size() == 1)
      return boost::any( 0 );
    else {
      std::ostringstream ss;
      ss << "No property '" << property << "' with " << pos.size() << " points implemented.";
      throw std::runtime_error(ss.str());
      return boost::any();
    }
  }

  /** @brief Getter for the lattice size */
  pos_t size() const { return L; }

  /** @brief Getter for the number of site types */
  int getMaxType() const override { return 1; }

private:

  /** @brief Prints the label of a given site */
  std::string site_label (int i) const {
      return "( " + boost::lexical_cast<std::string>(a * i) + " )";
  }

  /** @brief Prints the label of a given pair of sites */
  std::string bond_label (int i, int j) const {
    return (  "( " + boost::lexical_cast<std::string>(a * i) + " )"
            + " -- "
            + "( " + boost::lexical_cast<std::string>(a * j) + " )");
  }

private:
  int L;
  double a;
  bool pbc;
};

#endif
