/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef ADJACENCY_H
#define ADJACENCY_H

#include <vector>
#include <stdexcept>

namespace adj {
    class Adjacency
    {
    public:
        virtual ~Adjacency(){};
        virtual std::vector<int> forward(int) const = 0;
        virtual std::vector<int> all(int) const = 0;
        virtual int size() const = 0;
        
        virtual bool wraps_pbc(int, int) { return false; }
        
        virtual int site_type(int) { return 0; }
        virtual int bond_type(int, int) { return 0; }
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
            
            //        maquis::cout << p << " -> ";
            //        std::copy(ret.begin(), ret.end(), std::ostream_iterator<int>(maquis::cout, " "));
            //        maquis::cout << std::endl;
            
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
}

#endif
