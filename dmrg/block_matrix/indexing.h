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
#include <map>

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
public:
    typedef typename SymmGroup::charge charge;
    typedef typename std::vector<std::pair<typename SymmGroup::charge, std::size_t> >::iterator iterator;
    typedef typename std::vector<std::pair<typename SymmGroup::charge, std::size_t> >::const_iterator const_iterator;
    
    typedef basis_iterator_<SymmGroup> basis_iterator;
    
    std::size_t size_of_block(charge c) const // rename
    {
        assert( has(c) );
        return std::find_if(this->begin(), this->end(),
                            index_detail::is_first_equal<SymmGroup>(c))->second;
    }
    
    std::size_t position(charge c) const
    {
        return std::find_if(this->begin(), this->end(),
                            index_detail::is_first_equal<SymmGroup>(c)) - this->begin();
    }
    
    std::size_t destination(charge c) const
    {
        return std::find_if(this->begin(), this->end(),
                            boost::lambda::bind(index_detail::lt<SymmGroup>,
                                                boost::lambda::_1,
                                                std::make_pair(c, 0))) - this->begin();
    }
    
    bool has(charge c) const
    {
        return std::find_if(this->begin(), this->end(),
                            index_detail::is_first_equal<SymmGroup>(c)) != this->end();
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
    
    template<class Fusion>
    void init(Index<SymmGroup> const & a,
              Index<SymmGroup> const & b,
              Fusion f)
    {
        std::map<charge, size_t> block_begins;
        
        for (typename Index<SymmGroup>::const_iterator it1 = a.begin(); it1 != a.end(); ++it1)
            for (typename Index<SymmGroup>::const_iterator it2 = b.begin(); it2 != b.end(); ++it2)
            {
                charge pc = f(it1->first, it2->first);
                
                keys_.push_back(std::make_pair(it1->first, it2->first));
                vals_.push_back(block_begins[pc]);
                size_[pc] += it1->second * it2->second;
                block_begins[pc] += it1->second * it2->second;
            }
    }
    
    size_t operator()(charge a, charge b) const
    {
        assert( std::count(keys_.begin(), keys_.end(), std::make_pair(a, b)) > 0 );
        size_t pos = std::find(keys_.begin(), keys_.end(), std::make_pair(a, b))-keys_.begin();
        return vals_[pos];
    }
    
    size_t size(charge a, charge b) const
    {
        return size(a, b, static_cast<charge(*)(charge, charge)>(SymmGroup::fuse));
    }
    template<class Fusion>
    size_t size(charge a, charge b, Fusion f) const
    {
        charge pc = f(a, b);
        assert(size_.count(pc) > 0);
        return size_.find(pc)->second;
    }

    
    
private:
    std::map<charge, size_t> size_;
    std::vector<std::pair<charge, charge> > keys_;
    std::vector<size_t> vals_;
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
    std::transform(nc.begin(), nc.end(), nd.begin(), std::back_inserter(ret),
                   std::make_pair<charge, std::size_t>);
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
