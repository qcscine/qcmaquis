#ifndef ADJACENCY_H
#define ADJACENCY_H

class Adjacency
{
public:
	virtual ~Adjacency(){};
    virtual std::vector<std::size_t> forward(std::size_t) const = 0;
    virtual std::vector<std::size_t> all(std::size_t) const = 0;
    virtual std::size_t size() const = 0;
};

class ChainAdj : public Adjacency
{
public:
    ChainAdj(std::size_t L) : L_(L) { }
    
    std::vector<std::size_t> forward(std::size_t p) const
    {
        if (p < L_-1)
            return std::vector<std::size_t>(1, p+1);
        else
            return std::vector<std::size_t>();
    }
    
    std::vector<std::size_t> all(std::size_t p) const
    {
        std::vector<std::size_t> ret;
        if (p < L_-1)
            ret.push_back(p+1);
        if (p > 0)
            ret.push_back(p-1);
        return ret;
    }
    
    std::size_t size() const { return L_; }
    
private:
    std::size_t L_;
};

class PeriodicChainAdj : public Adjacency
{
public:
    PeriodicChainAdj(std::size_t L) : L_(L) { }
    
    std::vector<std::size_t> forward(std::size_t p) const
    {
        return std::vector<std::size_t>(1, (p+1) % L_);
    }
    
    std::vector<std::size_t> all(std::size_t p) const
    {
        std::vector<std::size_t> ret;
        ret.push_back((p+1)%L_);
        ret.push_back((p-1+L_)%L_);
        return ret;
    }
    
    std::size_t size() const { return L_; }
    
private:
    std::size_t L_;
};

class SquareAdj : public Adjacency
{
public:
    SquareAdj(std::size_t L, std::size_t W)
    : L_(L)
    , W_(W)
    { }
    
    /*
     0 4  8 12
     1 5  9 13
     2 6 10 14
     3 7 11 15
     */
    std::vector<std::size_t> forward(std::size_t p) const
    {
        std::vector<std::size_t> ret;
        if (p+1 < L_*W_ && (p+1) % W_ != 0)
            ret.push_back(p+1);
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);
        
        return ret;
    }
    
    std::vector<std::size_t> all(std::size_t p) const
    {
        std::vector<std::size_t> ret = forward(p);
        if (p >= 1 && p % W_ != 0)
            ret.push_back(p-1);
        if (p >= W_)
            ret.push_back(p-W_);
        
        return ret;
    }
    
    std::size_t size() const { return L_*W_; }
    
private:
    std::size_t L_, W_;
};

class CylinderAdj : public Adjacency
{
public:
    CylinderAdj(std::size_t L, std::size_t W)
    : L_(L)
    , W_(W)
    { }
    
    std::vector<std::size_t> forward(std::size_t p) const
    {
        std::vector<std::size_t> ret;
        if (p+1 < L_*W_ && (p+1) % W_ != 0)
            ret.push_back(p+1);
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);
        if (p+(W_-1) < L_*W_ && p % W_ == 0)
            ret.push_back(p+(W_-1));
        
        return ret;
    }
    
    std::vector<std::size_t> all(std::size_t) const
    {
        throw std::runtime_error("Not implemented.");
        return std::vector<std::size_t>();
    }
    
    std::size_t size() const { return L_*W_; }
    
private:
    std::size_t L_, W_;
};

#endif
