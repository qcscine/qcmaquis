/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Michele Dolfi <dolfim@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef CONTINUOUS_LATTICE_H
#define CONTINUOUS_LATTICE_H

#include "dmrg/models/lattice.h"

#include <sstream>
#include <boost/lexical_cast.hpp>

#include "dmrg/utils/BaseParameters.h"

class ContChain : public Lattice
{
public:
    typedef Lattice::pos_t pos_t;
    
    ContChain (BaseParameters & parms, bool pbc_=false)
    : L(parms.get<int>("L"))
    , N(parms.get<int>("Ndiscr"))
    , a(parms.get<double>("a"))
    , pbc(pbc_)
    {}
    
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L*N-1)
            ret.push_back(i+1);
        if (pbc && i == L*N-1)
            ret.push_back(0);
        return ret;
    }
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L*N-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        if (pbc && i == L*N-1)
            ret.push_back(0);
        if (pbc && i == 0)
            ret.push_back(L*N-1);
        return ret;
    }
    
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( site_label(pos[0]) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(pos[0], pos[1]) );
        else if (property == "type" && pos.size() == 1)
            return boost::any( pos[0]%N );
        else if (property == "type" && pos.size() == 2)
            return boost::any( 0 );
        else if (property == "x" && pos.size() == 1)
            return boost::any( a/N/2. + a/N * pos[0] );
        else if (property == "dx" && pos.size() == 1)
            return boost::any( a/N );
        else if (property == "dx" && pos.size() == 2)
            return boost::any( a/N * (pos[1]-pos[0]) );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L*N-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L*N-1 );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( (pos[0] < pos[1]) );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }
    
    pos_t size() const
    {
        return L*N;
    }
    
private:
    
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(a/N/2. + a/N * i) + " )";
    }
    
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(a/N/2. + a/N * i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(a/N/2. + a/N * j) + " )");
    }
    
private:
    int N, L;
    double a;
    bool pbc;
    
};

class MixedContChain : public Lattice
{
public:
    typedef Lattice::pos_t pos_t;
    
    MixedContChain (BaseParameters & parms1, int L1_, BaseParameters & parms2, int L2_, bool pbc_=false)
    : L1(L1_)
    , L2(L2_)
    , Lphys(parms1.get<double>("L"))
    , N1(parms1.get<int>("Ndiscr"))
    , a1(parms1.get<double>("a"))
    , N2(parms2.get<int>("Ndiscr"))
    , a2(parms2.get<double>("a"))
    , pbc(pbc_)
    {
        assert( parms1.get<double>("L") == parms1.get<double>("L") );
        if (pbc)
            throw std::runtime_error("Periodic boundary conditions are not implemented for the MixedChain.");
    }
    
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L1+L2-1)
            ret.push_back(i+1);
        if (pbc && i == L1+L2-1)
            ret.push_back(0);
        return ret;
    }
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L1+L2-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        if (pbc && i == L1+L2-1)
            ret.push_back(0);
        if (pbc && i == 0)
            ret.push_back(L1+L2-1);
        return ret;
    }
    
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
            return boost::any( get_x(pos[0]) );
        else if (property == "dx" && pos.size() == 1)
            return boost::any( get_dx(pos[0]) );
        else if (property == "dx" && pos.size() == 2)
            return boost::any( get_dx(pos[0], pos[1]) );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L1+L2-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L1+L2-1 );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( (pos[0] < pos[1]) );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }
    
    pos_t size() const
    {
        return L1+L2;
    }
    
private:
    
    double get_x (int i) const
    {
        int i1, i2;
        if (i < L1) {
            i1 = i;
            i2 = 0;
        } else {
            i1 = L1;
            i2 = i - L1;
        }
        
        return a1/N1 * i1 + a2/N2 * i2;
    }
    
    double get_dx (int i, int j) const
    {
        double dx;
        if (std::max(i, j) <= L1)
            dx = a1/N1 * (j-i);
        else
            dx = a2/N2 * (j-i);
        return dx;
    }
    double get_dx (int i) const
    {
        return (i<L1) ? a1/N1 : a2/N2;
    }

    
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>( get_x(i) ) + " )";
    }
    
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>( get_x(i) ) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>( get_x(j) ) + " )");
    }
    
private:
    int N1, N2, L1, L2;
    double a1, a2, Lphys;
    bool pbc;
    
};


class MixedContChain_c : public Lattice
{
public:
    typedef Lattice::pos_t pos_t;
    
    MixedContChain_c (BaseParameters & parms1, int L1_, BaseParameters & parms2, int L2_, bool pbc_=false)
    : L1(L1_)
    , L2(L2_)
    , Lphys(parms1.get<double>("L"))
    , N1(parms1.get<int>("Ndiscr"))
    , a1(parms1.get<double>("a"))
    , N2(parms2.get<int>("Ndiscr"))
    , a2(parms2.get<double>("a"))
    , pbc(pbc_)
    {
        assert( parms1.get<double>("L") == parms1.get<double>("L") );
        if (pbc)
            throw std::runtime_error("Periodic boundary conditions are not implemented for the MixedChain.");
    }
    
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L1+L2-1)
            ret.push_back(i+1);
        if (pbc && i == L1+L2-1)
            ret.push_back(0);
        return ret;
    }
    std::vector<pos_t> all(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L1+L2-1)
            ret.push_back(i+1);
        if (i > 0)
            ret.push_back(i-1);
        if (pbc && i == L1+L2-1)
            ret.push_back(0);
        if (pbc && i == 0)
            ret.push_back(L1+L2-1);
        return ret;
    }
    
    boost::any get_prop_(std::string const & property, std::vector<pos_t> const & pos) const
    {
        if (property == "label" && pos.size() == 1)
            return boost::any( site_label(pos[0]) );
        else if (property == "label" && pos.size() == 2)
            return boost::any( bond_label(pos[0], pos[1]) );
        //            else if (property == "type" && pos.size() == 1)
        //                return boost::any( 0 );
        //            else if (property == "type" && pos.size() == 2)
        //                return boost::any( 0 );
        else if (property == "x" && pos.size() == 1)
            return boost::any( get_x(pos[0]) );
        else if (property == "dx" && pos.size() == 1)
            return boost::any( get_dx(pos[0]) );
        else if (property == "dx" && pos.size() == 2)
            return boost::any( get_dx(pos[0], pos[1]) );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L1+L2-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L1+L2-1 );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( (pos[0] < pos[1]) );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }
    
    pos_t size() const
    {
        return L1+L2;
    }
    
private:
    
    double get_x (int i) const
    {
        int i1, i2;
        double middle = 0;
        double left = (L1 > 0) ? a1/N1/2. : a2/N2/2.;
        if (L1 == 0) {
            i1 = 0;
            i2 = i;
            middle = 0.;
        } else if (i < L1) {
            i1 = i;
            i2 = 0;
            middle = 0.;
        } else {
            i1 = L1-1;
            i2 = i - L1;
            middle = 3./2. * std::min(a1/N1, a2/N2);
        }
        
        return left + a1/N1 * i1 + middle + a2/N2 * i2;
    }
    double get_dx (int i) const
    {
        return (i<L1) ? a1/N1 : a2/N2;
    }
    
    double get_dx (int i, int j) const
    {
        if (L1 == 0)
            return a2/N2 * (j-i);
        else if (std::max(i, j) < L1 || L2 == 0)
            return a1/N1 * (j-i);
        else if (std::min(i,j) >= L1)
            return a2/N2 * (j-i);
        else
            return (j-i)/std::abs(j-i) * ( a1/N1 * (L1-1 - std::min(i,j))
                                          + 3./2. * std::min(a1/N1, a2/N2)
                                          + a2/N2 * (std::max(i,j) - L1) );
    }
    
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>( get_x(i) ) + " )";
    }
    
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>( get_x(i) ) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>( get_x(j) ) + " )");
    }
    
private:
    int N1, N2, L1, L2;
    double a1, a2, Lphys;
    bool pbc;
    
};

#endif
