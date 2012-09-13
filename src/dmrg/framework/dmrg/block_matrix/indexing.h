/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef TENSOR_INDEXING_H
#define TENSOR_INDEXING_H

#include <vector>
#include <algorithm>
#include <utility>

#include <boost/unordered_map.hpp>
#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>
#include <boost/lambda/bind.hpp>

#ifdef PYTHON_EXPORTS
#include <mp_tensors/wrappers.h>
#endif

#ifdef HAVE_ALPS_HDF5
#include <alps/hdf5.hpp>
#endif

namespace index_detail
{
    template<class SymmGroup>
    bool lt(std::pair<typename SymmGroup::charge, std::size_t> const & a,
            std::pair<typename SymmGroup::charge, std::size_t> const & b)
    {
        return a.first < b.first;
    }
    
    template<class SymmGroup>
    bool gt(std::pair<typename SymmGroup::charge, std::size_t> const & a,
            std::pair<typename SymmGroup::charge, std::size_t> const & b)
    {
        return a.first > b.first;
    }
    
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
: public std::vector<std::pair<typename SymmGroup::charge, std::size_t> >
{
    typedef std::vector<std::pair<typename SymmGroup::charge, std::size_t> > base;
public:
    typedef typename SymmGroup::charge charge;
    typedef typename std::vector<std::pair<typename SymmGroup::charge, std::size_t> >::iterator iterator;
    typedef typename std::vector<std::pair<typename SymmGroup::charge, std::size_t> >::const_iterator const_iterator;
    
    typedef basis_iterator_<SymmGroup> basis_iterator;
    
    std::size_t size_of_block(charge c) const // rename
    {
        assert( has(c) );
        return std::lower_bound(this->begin(), this->end(), std::make_pair(c,0), index_detail::gt<SymmGroup>)->second;
    }
    
    std::size_t position(charge c) const
    {
        const_iterator match = std::lower_bound(this->begin(), this->end(), std::make_pair(c,0), index_detail::gt<SymmGroup>);
        if(match != this->end() && (*match).first != c) match = this->end();
        return (match - this->begin());
    }

    std::size_t position(std::pair<charge, std::size_t> x) const
    {
        assert( has(x.first) );
        assert( x.second < size_of_block(x.first) );
        return x.second + std::accumulate(this->begin(),
                                          std::lower_bound(this->begin(), this->end(), x, index_detail::gt<SymmGroup>),
                                          0,
                                          boost::lambda::_1 + boost::lambda::bind(index_detail::get_second<SymmGroup>, boost::lambda::_2)
                                          );
    }

    std::size_t destination(charge c) const
    {
        return std::upper_bound(this->begin(), this->end(), std::make_pair(c,0), index_detail::gt<SymmGroup>) - this->begin();
    }
    
    bool has(charge c) const
    {
        return std::binary_search(this->begin(), this->end(), std::make_pair(c,0), index_detail::gt<SymmGroup>);
    }
    
    void sort()
    {
        std::sort(this->begin(), this->end(), index_detail::gt<SymmGroup>);
    }
    
    std::size_t insert(std::pair<charge, std::size_t> const & x)
    {
        std::size_t d = destination(x.first);
        std::vector<std::pair<charge, std::size_t> >::insert(this->begin() + d, x);
        return d;
    }
    
    void insert(std::size_t position, std::pair<charge, std::size_t> const & x)
    {
        std::vector<std::pair<charge, std::size_t> >::insert(this->begin() + position, x);
    }
    
    bool operator==(Index const & o) const
    {
        return (this->size() == o.size()) && std::equal(this->begin(), this->end(), o.begin());
    }
    
    basis_iterator basis_begin() const
    {
        assert( this->size() > 0 );
        return basis_iterator(*this);
    }
    
    std::vector<charge> charges() const
    {
        std::vector<charge> ret(this->size());
        for (std::size_t k = 0; k < this->size(); ++k) ret[k] = (*this)[k].first;
        return ret;
    }
    
    std::vector<std::size_t> sizes() const
    {
        std::vector<std::size_t> ret(this->size());
        for (std::size_t k = 0; k < this->size(); ++k) ret[k] = (*this)[k].second;
        return ret;
    }
    
    std::size_t sum_of_sizes() const
    {
        return std::accumulate(this->begin(), this->end(), 0,
                               boost::lambda::_1 + boost::lambda::bind(index_detail::get_second<SymmGroup>, boost::lambda::_2));
    }

private:
    void push_back(std::pair<charge, std::size_t> const & x){
        base::push_back(x);
    }

public:    
#ifdef PYTHON_EXPORTS
    std::size_t py_insert(wrapped_pair<SymmGroup> p)
    {
        return this->insert(p.data_);
    }
#endif /* PYTHON_EXPORTS */
    
#ifdef HAVE_ALPS_HDF5
    void load(alps::hdf5::archive & ar)
    {
        typedef std::vector<std::pair<typename SymmGroup::charge, std::size_t> > my_type;
        ar >> alps::make_pvp("Index",
                             static_cast<my_type&>(*this));
    }
    void save(alps::hdf5::archive & ar) const
    {
        typedef std::vector<std::pair<typename SymmGroup::charge, std::size_t> > my_type;
        ar << alps::make_pvp("Index",
                             static_cast<my_type const &>(*this));
    }
#endif
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
        //assert( std::count(keys_.begin(), keys_.end(), std::make_pair(a, b)) > 0 );
        size_t res = (*keys_vals_.find(std::make_pair(a,b))).second;
        return res;
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

#ifdef WIN32
template<class A, class B> std::pair<A, B> mypair(A & a, B & b) { return std::pair<A,B>(a,b); }
#endif

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
#ifndef WIN32
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
            if (ret.has(pdc))
                ret[ret.position(pdc)].second += ps;
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
