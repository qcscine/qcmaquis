#ifndef ADJACENCY_H
#define ADJACENCY_H

class Adjacency
{
public:
	virtual ~Adjacency(){};
    virtual std::vector<int> forward(int) const = 0;
    virtual std::vector<int> all(int) const = 0;
    virtual int size() const = 0;
};

class ChainAdj : public Adjacency
{
public:
    ChainAdj(int L) : L_(L) { }
    
    std::vector<int> forward(int p) const
    {
        if (p < L_-1)
            return std::vector<int>(1, p+1);
        else
            return std::vector<int>();
    }
    
    std::vector<int> all(int p) const
    {
        std::vector<int> ret;
        if (p < L_-1)
            ret.push_back(p+1);
        if (p > 0)
            ret.push_back(p-1);
        return ret;
    }
    
    int size() const { return L_; }
    
private:
    int L_;
};

class PeriodicChainAdj : public Adjacency
{
public:
    PeriodicChainAdj(int L) : L_(L) { }
    
    std::vector<int> forward(int p) const
    {
        return std::vector<int>(1, (p+1) % L_);
    }
    
    std::vector<int> all(int p) const
    {
        std::vector<int> ret;
        ret.push_back((p+1)%L_);
        ret.push_back((p-1+L_)%L_);
        return ret;
    }
    
    int size() const { return L_; }
    
private:
    int L_;
};

class PeriodicLadderAdj : public Adjacency
{
public:
    PeriodicLadderAdj(int L) : L_(L) { }
    
    std::vector<int> forward(int p) const
    {
        int N = 2*L_;
        
        std::vector<int> ret;
        if (p % 2 == 0)
            ret.push_back((p+1) % N);
        ret.push_back((p+2) % N);
        
        return ret;
    }
    
    std::vector<int> all(int p) const
    {
        int N = 2*L_;
        
        std::vector<int> ret;
        if (p % 2 == 0)
            ret.push_back((p+1) % N);
        else
            ret.push_back((p-1+N) % N);
        
        ret.push_back((p+2) % N);
        ret.push_back((p-2+N) % N);
        
        return ret;
    }
    
    int size() const { return 2*L_; }
    
private:
    int L_;
};

class SquareAdj : public Adjacency
{
public:
    SquareAdj(int L, int W)
    : L_(L)
    , W_(W)
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
        
//        std::cout << p << " -> ";
//        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(std::cout, " "));
//        std::cout << std::endl;
        
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
    
private:
    int L_, W_;
};

class SnakeSquareAdj : public Adjacency
{
public:
    SnakeSquareAdj(int L, int W) : L_(L), W_(W) { }
    
    /*
     3 4 11 12 19
     2 5 10 13 18
     1 6  9 14 17
     0 7  8 15 16 */
    
    std::vector<int> forward(int p) const
    {
        std::vector<int> ret;
        
#define wrap(pp) { if((pp) < L_*W_) ret.push_back(pp); }    
        wrap(p+1);
        
        if ((p+1) % W_ != 0)
        {
            int p_in_bump = p % (2*W_);
            if (p_in_bump / W_ == 0)
                wrap(p + 2*W_ - 1 - 2*p_in_bump)
            else
                wrap(p + 1 + 2*(2*W_-p_in_bump-1))
        }
#undef wrap
        
//        std::cout << p << " -> ";
//        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(std::cout, " "));
//        std::cout << std::endl;
        
        return ret;
    }
    
    std::vector<int> all(int p) const
    {
        throw std::runtime_error("Not implemented.");
    }
    
    int size() const { return L_*W_; }
    
private:
    int L_, W_;
};

class CylinderAdj : public Adjacency
{
public:
    CylinderAdj(int L, int W)
    : L_(L)
    , W_(W)
    { }
    
    std::vector<int> forward(int p) const
    {
        std::vector<int> ret;

        // down
        if ((p+1)%W_ != 0)
            ret.push_back(p+1);
        
        // right
        if (p+W_ < L_*W_)
            ret.push_back(p+W_);
        
        // around bad direction
        if ((p+1)%W_ == 0)
            ret.push_back(p-W_+1);
        
//        if (p+1 < L_*W_ && (p+1) % W_ != 0)
//            ret.push_back(p+1);
//        if (p+W_ < L_*W_)
//            ret.push_back(p+W_);
//        if (p+(W_-1) < L_*W_ && p % W_ == 0)
//            ret.push_back(p+(W_-1));
        
        return ret;
    }
    
    std::vector<int> all(int) const
    {
        throw std::runtime_error("Not implemented.");
        return std::vector<int>();
    }
    
    int size() const { return L_*W_; }
    
private:
    int L_, W_;
};

class PeriodicSquareLatticeAdj : public Adjacency
{
public:
    PeriodicSquareLatticeAdj(int L, int W)
    : L_(L)
    , W_(W)
    { }
    
    std::vector<int> forward(int p) const
    {
        int N = L_*W_;
        
        std::vector<int> ret;
        
        // down
        if ((p+1)%W_ != 0)
            ret.push_back(p+1);
        
        // right
        ret.push_back((p+W_)%N);
        
        // around bad direction
        if ((p+1)%W_ == 0)
            ret.push_back(p-W_+1);
        
        return ret;
    }
    
    std::vector<int> all(int) const
    {
        throw std::runtime_error("Not implemented.");
        return std::vector<int>();
    }
    
    int size() const { return L_*W_; }
    
private:
    int L_, W_;
};

#endif
