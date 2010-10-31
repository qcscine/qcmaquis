#ifndef TENSOR_INDEXING_H
#define TENSOR_INDEXING_H

#include <vector>
#include <algorithm>
#include <utility>

#include <boost/array.hpp>
#include <boost/lambda/lambda.hpp>

enum IndexName { alpha, sigma, beta, empty };

template<class T> boost::array<T,1> _(const T &a)
{ 
    boost::array<T,1> r = {a};
    return r;
}

inline boost::array<IndexName, 2> operator,(const IndexName a, const IndexName b)
{
	boost::array<IndexName, 2> ret;
	ret[0] = a;
	ret[1] = b;
	return ret;
}

template<class T, int in>
boost::array<T, in+1> operator,(const boost::array<T, in> a, const T b)
{
	boost::array<T, in+1> ret;
	for (int i = 0; i < in; i++)
		ret[i] = a[i];
	ret[in] = b;
	return ret;
}

template<class T, int in>
boost::array<T, in+1> operator,(const T a, const boost::array<T, in> b)
{
	boost::array<T, in+1> ret;
	ret[0] = a;
	for (int i = 0; i < in; i++)
		ret[i+1] = b[i];
	return ret;
}

template<class SymmGroup>
class Index : public std::vector<std::pair<typename SymmGroup::charge, std::size_t> >
{
public:
    typedef typename SymmGroup::charge charge;

    IndexName name;
    
    Index(IndexName n = empty)
    : std::vector<std::pair<charge, std::size_t> >()
    , name(n)
    { }
    
    Index(std::size_t M, IndexName n = empty)
    : std::vector<std::pair<charge, std::size_t> >(1, std::make_pair(SymmGroup::SingletCharge, M))
    , name(n)
    { }
    
    Index(std::vector<charge> const & charges, std::vector<std::size_t> const & sizes, IndexName n = empty)
    : name(n)
    {
        assert(sizes.size() == charges.size());
        std::transform(charges.begin(), charges.end(), sizes.begin(), std::back_inserter(*this),
            std::make_pair<charge, std::size_t>);
    }
    
    template<int L>
    Index(boost::array<charge, L> const & charges, boost::array<std::size_t, L> const & sizes, IndexName n = empty)
    : name(n)
    {
        std::transform(charges.begin(), charges.end(), sizes.begin(), std::back_inserter(*this),
            std::make_pair<charge, std::size_t>);
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

    bool operator==(const Index<SymmGroup> &o)
    {
        return name == o.name && this->size() == o.size() && std::equal(this->begin(), this->end(), o.begin());
    }
};

template<class SymmGroup>
Index<SymmGroup> transpose(Index<SymmGroup> const & inp)
{
    typedef typename SymmGroup::charge charge;
    
    std::vector<charge> oc = inp.charges(), nc = inp.charges();
    std::transform(nc.begin(), nc.end(), nc.begin(), std::negate<charge>());
    std::sort(nc.begin(), nc.end());
    
    std::vector<std::size_t> nd(inp.size()), od = inp.sizes();
    for (unsigned int i = 0; i < nd.size(); ++i)
        nd[i] = od[std::find(oc.begin(), oc.end(),
                                    -nc[i])-oc.begin()];
    
    return Index<SymmGroup>(nc, nd, inp.name);
}

template<class SymmGroup>
std::ostream& operator<<(std::ostream& os, Index<SymmGroup> const & idx)
{
    if (idx.name != empty)
        os << idx.name << std::endl;
    std::vector<typename SymmGroup::charge> charges = idx.charges();
    std::copy(charges.begin(), charges.end(), std::ostream_iterator<typename SymmGroup::charge>(os, " "));
    os << std::endl;
    
    std::vector<std::size_t> sizes = idx.sizes();
    std::copy(sizes.begin(), sizes.end(), std::ostream_iterator<std::size_t>(os, " "));
    os << std::endl;
    
    return os;
}

#endif
