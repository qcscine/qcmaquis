#ifndef MAQUIS_DMRG_lattice_hpp
#define MAQUIS_DMRG_lattice_hpp

#include "dmrg/models/lattice.h"

class ChainLattice : public Lattice
{
public:
    typedef Lattice::pos_t pos_t;
    
    ChainLattice (BaseParameters & parms, bool pbc_=false)
    : L(parms["L"])
    , a(parms["a"])
    , pbc(pbc_)
    { }

    ChainLattice (int L_, bool pbc_=false, double a_=1.)
    : L(L_)
    , a(a_)
    , pbc(pbc_)
    { }
    
    std::vector<pos_t> forward(pos_t i) const
    {
        std::vector<pos_t> ret;
        if (i < L-1)
            ret.push_back(i+1);
        if (pbc && i == L-1)
            ret.push_back(0);
        return ret;
    }
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
            return boost::any( a * pos[0] );
        else if (property == "at_open_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && (pos[0]==0 || pos[0]==L-1) );
        else if (property == "at_open_left_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==0 );
        else if (property == "at_open_right_boundary" && pos.size() == 1)
            return boost::any( (!pbc) && pos[0]==L-1 );
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
        return L;
    }
    
private:
    
    std::string site_label (int i) const
    {
        return "( " + boost::lexical_cast<std::string>(a * i) + " )";
    }
    
    std::string bond_label (int i, int j) const
    {
        return (  "( " + boost::lexical_cast<std::string>(a * i) + " )"
                + " -- "
                + "( " + boost::lexical_cast<std::string>(a * j) + " )");
    }
    
private:
    int L;
    double a;
    bool pbc;
    
};


class SquareLattice : public Lattice
{
public:
    SquareLattice(BaseParameters & parms)
    : L_(parms["L"])
    , W_(parms["W"])
    , a(parms["a"])
    { }
    
    /*
     0 4  8 12
     1 5  9 13
     2 6 10 14
     3 7 11 15
     */
    std::vector<int> forward(int p) const
    {
        std::vector<int> ret;
        if (p+1 < L_*W_ && (p+1) % W_ != 0)
            ret.push_back(p+1);
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);
        
        //        maquis::cout << p << " -> ";
        //        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(maquis::cout, " "));
        //        maquis::cout << std::endl;
        
        return ret;
    }
    
    std::vector<int> all(int p) const
    {
        std::vector<int> ret = forward(p);
        if (p >= 1 && p % W_ != 0)
            ret.push_back(p-1);
        if (p >= W_)
            ret.push_back(p-W_);
        
        return ret;
    }
    
    int size() const { return L_*W_; }
    
    
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
            return boost::any( x(pos[0]) );
        else if (property == "y" && pos.size() == 1)
            return boost::any( y(pos[0]) );
        else if (property == "wraps_pbc" && pos.size() == 2)
            return boost::any( false );
        else {
            std::ostringstream ss;
            ss << "No property '" << property << "' with " << pos.size() << " points implemented."; 
            throw std::runtime_error(ss.str());
            return boost::any();
        }
    }

    
    
private:
    
    double x (int i) const
    { return a * int(i/W_); }
    double y (int i) const
    { return a * (i%W_); }
    
    std::string site_label (int i) const
    {
        return "( " + ( boost::lexical_cast<std::string>(x(i))
                       + "," + boost::lexical_cast<std::string>(y(i)) ) + " )";
    }
    
    std::string bond_label (int i, int j) const
    {
        return (  "( " + ( boost::lexical_cast<std::string>(x(i))
                          + "," + boost::lexical_cast<std::string>(y(i)) ) + " )"
                + " -- "
                + "( " + ( boost::lexical_cast<std::string>(x(j))
                          + "," + boost::lexical_cast<std::string>(y(j)) ) + " )" );
    }        
    
    
    int L_, W_;
    double a;
};

#endif
