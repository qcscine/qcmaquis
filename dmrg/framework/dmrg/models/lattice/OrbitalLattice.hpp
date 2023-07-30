/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Michele Dolfi <dolfim@phys.ethz.ch>
 *                            Sebastian Keller <sebkelle@phys.ethz.ch>
 *               2020- by Robin Feldmann <robinfe@phys.chem.ethz.ch>
 *               2022- by Alberto Baiardi <abaiardi@ethz.ch>
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

#ifndef MAQUIS_DMRG_ORBITAL_LATTICE
#define MAQUIS_DMRG_ORBITAL_LATTICE

#include "dmrg/models/lattice/lattice.h"
#include <sstream>
#include <vector>
#include <set>
#include <boost/lexical_cast.hpp>
#include <boost/lambda/lambda.hpp>
#include <numeric>
#include "dmrg/utils/BaseParameters.h"
#include "LatticeHelperClass.hpp"

class Orbitals : public lattice_impl
{
public:
  typedef int subcharge;
  typedef lattice_impl::pos_t pos_t;

  /** @brief Class contructor */
  Orbitals (BaseParameters & parms) : L(parms["L"]), irreps(L, 0), order(L)
  {
    // Extracts the order
    order = LatticeHelperClass::getOrbitalOrder(parms, "orbital_order", true);
    if (parms.is_set("integral_file")) {
      std::string integral_file = parms["integral_file"];
      if (!boost::filesystem::exists(integral_file))
        throw std::runtime_error("integral_file " + integral_file + " does not exist\n");
      std::ifstream orb_file;
      orb_file.open(parms["integral_file"].c_str());
      orb_file.ignore(std::numeric_limits<std::streamsize>::max(),'\n');
      std::string line;
      std::getline(orb_file, line);
      orb_file.close();
      std::vector<std::string> split_line;
      boost::split(split_line, line, boost::is_any_of("="));
      // record the site_types in parameters
      parms.set("site_types", split_line[1]);
      irreps = parse_irreps(split_line[1]);
    }
    else if (parms.is_set("site_types")) {
      std::vector<subcharge> symm_vec = parms["site_types"];
      assert(L == symm_vec.size());
      for (subcharge p = 0; p < L; ++p)
        irreps[p] = symm_vec[order[p]];
    }
    else
      throw std::runtime_error("\"integral_file\" in parms input file or site_types is not set\n");
    numTypes = *std::max_element(irreps.begin(), irreps.end())+1;
  }

  std::vector<pos_t> forward(pos_t i) const override {
      std::vector<pos_t> ret;
      if (i < L-1)
          ret.push_back(i+1);
      return ret;
  }

  std::vector<pos_t> all(pos_t i) const override {
      std::vector<pos_t> ret;
      if (i < L-1)
          ret.push_back(i+1);
      if (i > 0)
          ret.push_back(i-1);
      return ret;
  }

  boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const override
  {
    if (property == "label" && pos.size() == 1) // return "( label )" as string
      return boost::any( site_label(order[pos[0]]) );
    if (property == "label_int" && pos.size() == 1) // just return label as integer
      return boost::any(order[pos[0]]);
    else if (property == "label" && pos.size() == 2)
      return boost::any( bond_label(order[pos[0]], order[pos[1]]) );
    else if (property == "type" && pos.size() == 1)
      return boost::any( irreps[pos[0]] );
    else if (property == "type" && pos.size() == 2)
      return boost::any( 0 );
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
  pos_t size() const override { return L; }

  /** @brief Getter for the number of types of sites */
  int getMaxType() const override { return numTypes; }

private:

    std::vector<subcharge> parse_irreps(std::string input)
    {
        std::vector<subcharge> symm_vec, ret(L, 0);

        std::replace(input.begin(), input.end(), ',', ' ');
        std::istringstream iss(input);
        subcharge number;
        while( iss >> number )
            symm_vec.push_back(number-1);

        assert(L == symm_vec.size());
        for (subcharge p = 0; p < L; ++p)
            ret[p] = symm_vec[order[p]];

        return ret;
    }

    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(i) + " )";
    }

    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(j) + " )");
    }

private:
    pos_t L;
    int numTypes;
    std::vector<subcharge> irreps;
    std::vector<int> order;
};

#endif
