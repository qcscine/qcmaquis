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

#ifndef TENSOR_INDEXING_H
#define TENSOR_INDEXING_H

#include <vector>
#include <algorithm>
#include <utility>

#include <boost/unordered_map.hpp>
#include <boost/container/flat_set.hpp>
#include <boost/container/flat_map.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#ifdef PYTHON_EXPORTS
#include <mp_tensors/wrappers.h>
#endif

#include <boost/serialization/serialization.hpp>
#include <boost/serialization/utility.hpp>
#include <boost/serialization/vector.hpp>


namespace index_detail
{
    template<class SymmGroup>
    bool lt(std::pair<typename SymmGroup::charge, std::size_t> const & a,
        std::pair<typename SymmGroup::charge, std::size_t> const & b)
    {
        return a.first < b.first;
    }  

    template<class SymmGroup>
    struct gt{
        bool operator()(std::pair<typename SymmGroup::charge, std::size_t> const & a,
                        std::pair<typename SymmGroup::charge, std::size_t> const & b){
            return a.first > b.first;
        }
    };

    template<class SymmGroup>
    typename SymmGroup::charge get_first(std::pair<typename SymmGroup::charge, std::size_t> const & x)
    {
        return x.first;
    }
    
    template<class SymmGroup>
    std::size_t get_second(std::pair<typename SymmGroup::charge, std::size_t> const & x)
    {
        return x.second;
    }
    
    // simpler, and potentially faster since inlining is easier for the compiler
    template<class SymmGroup>
    class is_first_equal
    {
    public:
        is_first_equal(typename SymmGroup::charge c) : c_(c) { }
        
        bool operator()(std::pair<typename SymmGroup::charge, std::size_t> const & x) const
        {
            return x.first == c_;
        }
        
    private:
        typename SymmGroup::charge c_;
    };
}

template<class SymmGroup>
class basis_iterator_;

template<class SymmGroup> class Index
{
    typedef std::vector<std::pair<typename SymmGroup::charge, std::size_t> > data_type;
    
public:
    typedef typename SymmGroup::charge charge;
    typedef typename data_type::value_type value_type;
    
    typedef typename data_type::iterator iterator;
    typedef typename data_type::const_iterator const_iterator;
    
    typedef typename data_type::reverse_iterator reverse_iterator;
    typedef typename data_type::const_reverse_iterator const_reverse_iterator;
    
    typedef basis_iterator_<SymmGroup> basis_iterator;
    
    Index() : sorted_(true) {}
    
    std::size_t size_of_block(charge c) const
    {
        assert( has(c) );
        return (*this)[position(c)].second;
    }

    std::size_t size_of_block(charge c, bool position_check) const
    {
        // I have to ignore the position_check argument because I can't dereference the end() iterator anyway
        std::size_t pos = position(c);
        if (pos == data_.size())
            return 0;
        return (*this)[pos].second;
    }
    
    std::size_t position(charge c) const
    {
        const_iterator match;
        if (sorted_)
            match = std::lower_bound(data_.begin(), data_.end(), std::make_pair(c,0), index_detail::gt<SymmGroup>());
        else
            match = std::find_if(data_.begin(), data_.end(), index_detail::is_first_equal<SymmGroup>(c));
        
        if (match != data_.end() && (*match).first != c) match = data_.end();
        return std::distance(data_.begin(), match);
    }

    std::size_t position(std::pair<charge, std::size_t> x) const
    {
        assert( has(x.first) );
        assert( x.second < size_of_block(x.first) );
        const_iterator to = data_.begin()+position(x.first);
        return x.second + std::accumulate(data_.begin(), to, 0,
                                          boost::lambda::_1 + boost::lambda::bind(index_detail::get_second<SymmGroup>, boost::lambda::_2)
                                         );
    }
    
    bool has(charge c) const
    {
        if (sorted_)
            return std::binary_search(data_.begin(), data_.end(), std::make_pair(c,0), index_detail::gt<SymmGroup>());
        else
            return std::find_if(data_.begin(), data_.end(),
                                index_detail::is_first_equal<SymmGroup>(c)) != data_.end();
    }
    
    void sort()
    {
        std::sort(data_.begin(), data_.end(), index_detail::gt<SymmGroup>());
        sorted_ = true;
    }
    
    std::size_t insert(std::pair<charge, std::size_t> const & x)
    {
        if (sorted_) {
            std::size_t d = destination(x.first);
            data_.insert(data_.begin() + d, x);
            return d;
        } else {
            push_back(x);
            return data_.size()-1;
        }
    }
    
    void insert(std::size_t position, std::pair<charge, std::size_t> const & x)
    {
        data_.insert(data_.begin() + position, x);
        sorted_ = false;
    }
    
    void shift(charge diff)
    {
        for (std::size_t k = 0; k < data_.size(); ++k)
            (*this)[k].first = SymmGroup::fuse((*this)[k].first, diff);
    }
    
    bool operator==(Index const & o) const
    {
        return (data_.size() == o.size()) && std::equal(data_.begin(), data_.end(), o.begin());
    }

    bool operator!=(Index const & o) const
    {
        return !( *this == o );
    }

    basis_iterator basis_begin() const
    {
        assert( data_.size() > 0 );
        return basis_iterator(*this);
    }
    
    std::vector<charge> charges() const
    {
        std::vector<charge> ret(data_.size());
        for (std::size_t k = 0; k < data_.size(); ++k) ret[k] = (*this)[k].first;
        return ret;
    }
    
    std::vector<std::size_t> sizes() const
    {
        std::vector<std::size_t> ret(data_.size());
        for (std::size_t k = 0; k < data_.size(); ++k) ret[k] = (*this)[k].second;
        return ret;
    }
    
    std::size_t sum_of_sizes() const
    {
		//boost::function<std::size_t (std::size_t,std::size_t)> pred = boost::lambda::_1 + boost::lambda::bind(index_detail::get_second<SymmGroup>, boost::lambda::_2);
        return std::accumulate(data_.begin(), data_.end(), 0, boost::lambda::ret<std::size_t>(boost::lambda::_1 + boost::lambda::bind(index_detail::get_second<SymmGroup>, boost::lambda::_2)));
    }

    // This is mostly forwarding of the std::vector
    iterator begin() { return data_.begin(); }
    iterator end() { return data_.end(); }
    const_iterator begin() const { return data_.begin(); }
    const_iterator end() const { return data_.end(); }
    
    reverse_iterator rbegin() { return data_.rbegin(); }
    reverse_iterator rend() { return data_.rend(); }
    const_reverse_iterator rbegin() const { return data_.rbegin(); }
    const_reverse_iterator rend() const { return data_.rend(); }
    
    value_type & operator[](std::size_t p) { return data_[p]; }
    value_type const & operator[](std::size_t p) const { return data_[p]; }
    
    boost::tuple<charge, std::size_t> element(std::size_t p) const
    {
        std::size_t i=0;
        while (p >= (*this)[i].second) {
            p -= (*this)[i].second;
            ++i;
        }
        return boost::make_tuple( (*this)[i].first, p );
    }

    std::size_t size() const { return data_.size(); }
    
    iterator erase(iterator p) { iterator r = data_.erase(p); return r; }
    iterator erase(iterator a, iterator b) { iterator r = data_.erase(a,b); return r; }
    
    friend void swap(Index & a, Index & b)
    {
        using std::swap;
        swap(a.data_,   b.data_);
        swap(a.sorted_, b.sorted_);
    }
    
private:
    data_type data_;
    bool sorted_;
    
    void push_back(std::pair<charge, std::size_t> const & x){
        data_.push_back(x);
    }
    
    std::size_t destination(charge c) const
    {
        return std::find_if(data_.begin(), data_.end(),
                            boost::lambda::bind(index_detail::lt<SymmGroup>,
                                                boost::lambda::_1,
                                                std::make_pair(c, 0))) - data_.begin();
    }

public:
#ifdef PYTHON_EXPORTS
    std::size_t py_insert(wrapped_pair<SymmGroup> p)
    {
        return data_.insert(p.data_);
    }
#endif /* PYTHON_EXPORTS */
   
    template <class Archive>
    void load(Archive & ar)
    {
        ar["Index"] >> data_;
    }
    template <class Archive>
    void save(Archive & ar) const
    {
        ar["Index"] << data_;
    }
    
    friend class boost::serialization::access;

    template <class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar & data_;
    }
    template <class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & data_;
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template<class SymmGroup>
class ProductBasis
{
public:
    typedef typename SymmGroup::charge charge;
    typedef std::size_t size_t;
    
    ProductBasis(Index<SymmGroup> const & a,
                 Index<SymmGroup> const & b)
    {
        init(a, b, static_cast<charge(*)(charge, charge)>(SymmGroup::fuse));
    }
    
    template<class Fusion>
    ProductBasis(Index<SymmGroup> const & a,
                 Index<SymmGroup> const & b,
                 Fusion f)
    {
        init(a, b, f);
    }
   
private: 
    template<class Fusion>
    void init(Index<SymmGroup> const & a,
              Index<SymmGroup> const & b,
              Fusion f)
    {
        keys_vals_.rehash((keys_vals_.size() + a.size()*b.size()) / keys_vals_.max_load_factor() + 1); // from http://www.boost.org/doc/libs/1_37_0/doc/html/unordered/buckets.html
        for (typename Index<SymmGroup>::const_iterator it1 = a.begin(); it1 != a.end(); ++it1)
            for (typename Index<SymmGroup>::const_iterator it2 = b.begin(); it2 != b.end(); ++it2)
            {
                charge pc = f(it1->first, it2->first);
                keys_vals_[std::make_pair(it1->first, it2->first)] = size_[pc]; 	
          //    keys_vals_.insert(std::make_pair(std::make_pair(it1->first, it2->first),size_[pc])); 	
                size_[pc] += it1->second * it2->second;
            }
    }

public:
    size_t operator()(charge a, charge b) const
    {
        typedef typename boost::unordered_map<std::pair<charge, charge>, size_t>::const_iterator match_type;
        match_type match = keys_vals_.find(std::make_pair(a,b));
        assert( match != keys_vals_.end() );
        return match->second;
    }
    
    inline size_t size(charge pc) const
    {
        assert(size_.count(pc) > 0);
        return size_[pc];
    }
    
    // for the moment let's avoid the default template argument (C++11)
    inline size_t size(charge a, charge b) const
    {
        return size(a, b, static_cast<charge(*)(charge, charge)>(SymmGroup::fuse));
    }
    template<class Fusion>
    size_t size(charge a, charge b, Fusion f) const
    {
        charge pc = f(a, b);
        assert(size_.count(pc) > 0);
        return size_[pc];
    }

private:
    mutable boost::unordered_map<charge, size_t> size_;
    boost::unordered_map<std::pair<charge, charge>, size_t> keys_vals_;
};

template<class SymmGroup>
class basis_iterator_
{
public:
    typedef typename SymmGroup::charge charge;
    
    basis_iterator_(Index<SymmGroup> const & idx, bool at_end = false)
    : idx_(idx)
    , cur_block(idx.begin())
    , cur_i(0)
    , max_i(cur_block->second)
    { }
    
    std::pair<charge, std::size_t> operator*() const
    {
        return std::make_pair(cur_block->first, cur_i);
    }
    
    boost::shared_ptr<std::pair<charge, std::size_t> > operator->() const
    {
        return boost::shared_ptr<std::pair<charge, std::size_t> >(new std::pair<charge, std::size_t>(cur_block->first, cur_i));
    }
    
    basis_iterator_ & operator++()
    {
        ++cur_i;
        if (cur_i != max_i)
            return *this;
        else {
            ++cur_block;
            if (cur_block != idx_.end()) {
                cur_i = 0;
                max_i = cur_block->second;
            }
            return *this;
        }
    }
    
    basis_iterator_ operator+(int k)
    {
        assert( k >= 0 );
        basis_iterator_ r = *this;
        for ( ; k > 0; --k)
            ++r;
        return r;
    }
    
    bool end() const
    {
        return cur_block == idx_.end();
    }
    
private:
    Index<SymmGroup> const & idx_;
    typename Index<SymmGroup>::const_iterator cur_block;
    std::size_t cur_i, max_i;
};

template<class SymmGroup>
basis_iterator_<SymmGroup> operator+(basis_iterator_<SymmGroup> it, std::size_t p)
{
    for ( ; p > 0; --p)
        ++it;
    return it;
}

// This is a workaround for MSVC++
// It can be removed as soon as this bug is fixed:
// http://social.msdn.microsoft.com/Forums/en/vclanguage/thread/bab04536-8a8d-4b5e-9a49-e10144688667

#if defined(WIN32) || defined(WIN64)
template<class A, class B> std::pair<A, B> mypair(A & a, B & b) { return std::pair<A,B>(a,b); }
#endif

// with sorted index we actually impose strong equality
template<class SymmGroup>
bool weak_equal(Index<SymmGroup> const & a, Index<SymmGroup> const & b)
{
    return (a == b);
}

template<class SymmGroup>
Index<SymmGroup> adjoin(Index<SymmGroup> const & inp)
{
    typedef typename SymmGroup::charge charge;
    
    std::vector<charge> oc = inp.charges(), nc = inp.charges();
    std::transform(nc.begin(), nc.end(), nc.begin(), std::negate<charge>());
    std::sort(nc.begin(), nc.end());
    
    std::vector<std::size_t> nd(inp.size()), od = inp.sizes();
    for (unsigned int i = 0; i < nd.size(); ++i)
        nd[i] = od[std::find(oc.begin(), oc.end(),
                             -nc[i])-oc.begin()];
    
    Index<SymmGroup> ret;
    for (size_t i=0; i<nc.size(); ++i)
#if not defined(WIN32) && not defined(WIN64)
        ret.insert(std::make_pair(nc[i], nd[i]));
#else
        ret.insert(mypair(nc[i], nd[i]));
#endif
    return ret;
}   

template<class SymmGroup>
std::ostream& operator<<(std::ostream& os, Index<SymmGroup> const & idx)
{
    os << "|";
    for (typename Index<SymmGroup>::const_iterator it = idx.begin();
         it != idx.end();
         ++it)
    {
        os << "( " << it->first << ": " << it->second << " )";
    }
    os << "|";
    
    return os;
}

template<class SymmGroup>
Index<SymmGroup> operator*(Index<SymmGroup> const & i1,
                           Index<SymmGroup> const & i2)
{
    typedef typename SymmGroup::charge charge;
    
    Index<SymmGroup> ret;
    for (typename Index<SymmGroup>::const_iterator it1 = i1.begin(); it1 != i1.end(); ++it1)
        for (typename Index<SymmGroup>::const_iterator it2 = i2.begin(); it2 != i2.end(); ++it2)
        {
            charge pdc = SymmGroup::fuse(it1->first, it2->first);
            std::size_t ps = it1->second * it2->second;
            std::size_t match = ret.position(pdc);
            if (match < ret.size())
                ret[match].second += ps;
            else
                ret.insert(std::make_pair(pdc, ps));
        }
    ret.sort();
    return ret;
}

template<class SymmGroup>
Index<SymmGroup> common_subset(Index<SymmGroup> & a,
                               Index<SymmGroup> & b)
{
    a.erase(std::remove_if(a.begin(), a.end(),
                           !boost::lambda::bind(&Index<SymmGroup>::has, b,
                                                boost::lambda::bind(index_detail::get_first<SymmGroup>, boost::lambda::_1))),
            a.end());
    
    b.erase(std::remove_if(b.begin(), b.end(),
                           !boost::lambda::bind(&Index<SymmGroup>::has, a,
                                                boost::lambda::bind(index_detail::get_first<SymmGroup>, boost::lambda::_1))),
            b.end());
    return a;
}

template<class charge>
std::pair<charge, std::size_t> operator-(std::pair<charge, std::size_t> const & p)
{
    return std::make_pair(-p.first, p.second);
}

template<class T> boost::array<T, 1> _(T const & a)
{ 
    boost::array<T, 1> r;
    r[0] = a;
    return r;
}

template<class T> boost::array<T, 2> _(T const & a, T const & b)
{ 
    boost::array<T, 2> r;
    r[0] = a;
    r[1] = b;
    return r;
}

#define IMPL_COMMA(tpl, type) \
tpl boost::array<type, 2> operator^(type const & a, type const & b) { \
    boost::array<type, 2> ret; \
    ret[0] = a; \
    ret[1] = b; \
    return ret; \
}
#define IMPL_COMMA_2(tpl, type) \
tpl boost::array<type, L+1> operator^(boost::array<type, L> const & a, type const & b) { \
    boost::array<type, L+1> ret; \
    std::copy(a.begin(), a.end(), ret.begin()); \
    ret[L] = b; \
    return ret; \
}

#define CO ,

IMPL_COMMA(template<class SymmGroup>, Index<SymmGroup>)
IMPL_COMMA(template<class charge>, std::pair<charge CO std::size_t>)

#undef CO
#undef IMPL_COMMA
#undef IMPL_COMMA_2

template<class T, unsigned long L>
boost::array<T, L+1> operator^(boost::array<T, L> const & a, T const & b)
{
	boost::array<T, L+1> ret;
    std::copy(a.begin(), a.end(), ret.begin());
	ret[L] = b;
	return ret;
}

template<class T, unsigned long L>
boost::array<T, L+1> operator^(T const & a, boost::array<T, L> const & b)
{
	boost::array<T, L+1> ret;
	ret[0] = a;
	for (int i = 0; i < L; i++)
		ret[i+1] = b[i];
	return ret;
}


#endif
