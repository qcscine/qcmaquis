/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *               2014-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef TENSOR_DUAL_INDEX_H
#define TENSOR_DUAL_INDEX_H

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/serialization/nvp.hpp>
#include <boost/preprocessor/repetition.hpp>

namespace dual_index_detail
{
    template <class SymmGroup>
    struct QnBlock
    {
        typename SymmGroup::charge lc;
        typename SymmGroup::charge rc;
        std::size_t                ls;
        std::size_t                rs;

        bool operator==(QnBlock const & o) const
        {
            return lc == o.lc && rc == o.rc && ls == o.ls && rs == o.rs;
        }
    };

    template <class SymmGroup>
    QnBlock<SymmGroup> make_QnBlock(typename SymmGroup::charge lc, typename SymmGroup::charge rc,
                                    std::size_t ls, std::size_t rs)
    {
        QnBlock<SymmGroup> ret;
        ret.lc = lc;
        ret.rc = rc;
        ret.ls = ls;
        ret.rs = rs;
        return ret;
    }

    template<class SymmGroup>
    struct gt {
        bool operator()(QnBlock<SymmGroup> const & a,
                        QnBlock<SymmGroup> const & b)
        {
            if (a.lc > b.lc)
                return true;
            else if (a.lc < b.lc)
                return false;
            else
                return a.rc > b.rc;
        }
    };

    template<class SymmGroup>
    struct gt_row{
        bool operator()(QnBlock<SymmGroup> const & a,
                        QnBlock<SymmGroup> const & b)
        {
            return (a.lc > b.lc);
        }
    };

    template<class SymmGroup>
    bool lt(QnBlock<SymmGroup> const & a,
            QnBlock<SymmGroup> const & b)
    {
        if (a.lc < b.lc)
            return true;
        else if (a.lc > b.lc)
            return false;
        else
            return a.rc < b.rc;
    }

    //// simpler, and potentially faster since inlining is easier for the compiler
    template<class SymmGroup>
    class is_first_equal
    {
    public:
        is_first_equal(typename SymmGroup::charge c1, typename SymmGroup::charge c2) : c1_(c1), c2_(c2) { }

        bool operator()(QnBlock<SymmGroup> const & x) const
        {
            return x.lc == c1_ && x.rc == c2_;
        }

    private:
        typename SymmGroup::charge c1_;
        typename SymmGroup::charge c2_;
    };
}

/*
namespace boost { namespace serialization {
// from http://uint32t.blogspot.ch/2008/03/update-serializing-boosttuple-using.html

#define GENERATE_ELEMENT_SERIALIZE(z,which,unused) \
    ar & boost::serialization::make_nvp("element",t.get< which >());

#define GENERATE_TUPLE_SERIALIZE(z,nargs,unused)                        \
    template< typename Archive, BOOST_PP_ENUM_PARAMS(nargs,typename T) > \
    void serialize(Archive & ar,                                        \
                   boost::tuple< BOOST_PP_ENUM_PARAMS(nargs,T) > & t,   \
                   const unsigned int version)                          \
    {                                                                   \
      BOOST_PP_REPEAT_FROM_TO(0,nargs,GENERATE_ELEMENT_SERIALIZE,~);    \
    }

    BOOST_PP_REPEAT_FROM_TO(1,6,GENERATE_TUPLE_SERIALIZE,~);
}}
*/

namespace boost { namespace serialization {

    template <class Archive, class SymmGroup>
    void serialize(Archive & ar, dual_index_detail::QnBlock<SymmGroup> & t,
                   const unsigned int version)
    {
        ar & boost::serialization::make_nvp("element",t.lc);
        ar & boost::serialization::make_nvp("element",t.rc);
        ar & boost::serialization::make_nvp("element",t.ls);
        ar & boost::serialization::make_nvp("element",t.rs);
    }

}}


template<class SymmGroup> class DualIndex
: protected std::vector<dual_index_detail::QnBlock<SymmGroup> >
{
    typedef std::vector<dual_index_detail::QnBlock<SymmGroup> > base_t;
    
public:
    typedef typename SymmGroup::charge charge;
    typedef typename base_t::value_type value_type;
    
    typedef typename base_t::iterator iterator;
    typedef typename base_t::const_iterator const_iterator;
    
    typedef typename base_t::reverse_iterator reverse_iterator;
    typedef typename base_t::const_reverse_iterator const_reverse_iterator;
    
    typedef basis_iterator_<SymmGroup> basis_iterator;
    
    DualIndex() : sorted_(true) {}
    
    //std::size_t size_of_block(charge c) const
    //{
    //    assert( has(c) );
    //    return (*this)[position(c)].second;
    //}

    //std::size_t size_of_block(charge c, bool position_check) const
    //{
    //    // I have to ignore the position_check argument because I can't dereference the end() iterator anyway
    //    std::size_t pos = position(c);
    //    if (pos == this->size())
    //        return 0;
    //    return (*this)[pos].second;
    //}

    std::size_t left_block_size(charge r, charge c) const {
        std::size_t pos = position(dual_index_detail::make_QnBlock<SymmGroup>(r,c,0,0));
        return (*this)[pos].ls;
    }
    std::size_t right_block_size(charge r, charge c) const {
        std::size_t pos = position(dual_index_detail::make_QnBlock<SymmGroup>(r,c,0,0));
        return (*this)[pos].rs;
    }
    
    std::size_t position(charge row, charge col) const
    {
        const_iterator match;
        if (sorted_)
            match = std::lower_bound(this->begin(), this->end(), dual_index_detail::make_QnBlock<SymmGroup>(row,col,std::size_t(0),std::size_t(0)), dual_index_detail::gt<SymmGroup>());
        else
            match = std::find_if(this->begin(), this->end(), dual_index_detail::is_first_equal<SymmGroup>(row,col));
        
        if (match != this->end() && ((*match).lc != row || (*match).rc != col))
            match = this->end();
        return std::distance(this->begin(), match);
    }

    bool has(charge row, charge col) const
    {
        if (sorted_)
            return std::binary_search(this->begin(), this->end(), dual_index_detail::make_QnBlock<SymmGroup>(row,col,0,0), dual_index_detail::gt<SymmGroup>());
        else
            return std::find_if(this->begin(), this->end(),
                                dual_index_detail::is_first_equal<SymmGroup>(row,col)) != this->end();
    }

    bool is_sorted() const { return sorted_; }
    
    void sort()
    {
        std::sort(this->begin(), this->end(), dual_index_detail::gt<SymmGroup>());
        sorted_ = true;
    }
    
    std::size_t insert(value_type const & x)
    {
        if (sorted_) {
            std::size_t d = destination(x);
            base_t::insert(this->begin() + d, x);
            return d;
        } else {
            push_back(x);
            return this->size()-1;
        }
    }
    
    void insert(std::size_t position, std::pair<charge, std::size_t> const & x)
    {
        base_t::insert(this->begin() + position, x);
        sorted_ = false;
    }
    
    void shift(charge diff)
    {
        for (std::size_t k = 0; k < this->size(); ++k)
        {
            (*this)[k].lc = SymmGroup::fuse((*this)[k].lc, diff);
            (*this)[k].rc = SymmGroup::fuse((*this)[k].rc, diff);
        }
    }
    
    bool operator==(DualIndex const & o) const
    {
        return (this->size() == o.size()) && std::equal(this->begin(), this->end(), o.begin());
    }

    bool operator!=(DualIndex const & o) const
    {
        return !( *this == o );
    }

    basis_iterator basis_begin() const
    {
        assert( this->size() > 0 );
        return basis_iterator(*this);
    }
    
    /*
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
    */
    
    std::size_t sum_of_left_sizes() const
    {
        return std::accumulate(this->begin(), this->end(), 0,
                               boost::lambda::_1
                               //+ boost::lambda::bind(static_cast<std::size_t(*)(value_type)>(boost::tuples::get<2>), boost::lambda::_2)
                               + boost::lambda::bind(&dual_index_detail::QnBlock<SymmGroup>::ls, boost::lambda::_2)
                              );
    }

    // This is mostly forwarding of the std::vector
    iterator begin() { return base_t::begin(); }
    iterator end() { return base_t::end(); }
    const_iterator begin() const { return base_t::begin(); }
    const_iterator end() const { return base_t::end(); }
    
    reverse_iterator rbegin() { return base_t::rbegin(); }
    reverse_iterator rend() { return base_t::rend(); }
    const_reverse_iterator rbegin() const { return base_t::rbegin(); }
    const_reverse_iterator rend() const { return base_t::rend(); }

    void resize(std::size_t sz) { base_t::resize(sz); }
    
    value_type & operator[](std::size_t p) { return static_cast<base_t&>(*this)[p]; }
    value_type const & operator[](std::size_t p) const { return static_cast<base_t const&>(*this)[p]; }
    
    //boost::tuple<charge, std::size_t> element(std::size_t p) const
    //{
    //    std::size_t i=0;
    //    while (p >= (*this)[i].second) {
    //        p -= (*this)[i].second;
    //        ++i;
    //    }
    //    return dual_index_detail::make_QnBlock<SymmGroup>( (*this)[i].first, p );
    //}

    std::size_t size() const { return base_t::size(); }
    
    iterator erase(iterator p) { iterator r = base_t::erase(p); return r; }
    iterator erase(iterator a, iterator b) { iterator r = base_t::erase(a,b); return r; }
    
private:
    bool sorted_;
    
    void push_back(value_type const & x){
        base_t::push_back(x);
    }
    
    std::size_t destination(value_type const & x) const
    {
        return std::find_if(this->begin(), this->end(),
                            boost::lambda::bind(dual_index_detail::lt<SymmGroup>,
                                                boost::lambda::_1,
                                                x)) - this->begin();
    }

public:
#ifdef PYTHON_EXPORTS
    std::size_t py_insert(wrapped_pair<SymmGroup> p)
    {
        return this->insert(p.data_);
    }
#endif /* PYTHON_EXPORTS */
   
    template <class Archive>
    void load(Archive & ar)
    {
        typedef std::vector<boost::tuple<typename SymmGroup::charge, typename SymmGroup::charge, std::size_t, std::size_t> > my_type;
        ar["DualIndex"] >> static_cast<my_type&>(*this);
    }
    template <class Archive>
    void save(Archive & ar) const
    {
        typedef std::vector<std::pair<typename SymmGroup::charge, std::size_t> > my_type;
        ar["DualIndex"] << static_cast<value_type const &>(*this);
    }
    
    friend class boost::serialization::access;

    template <class Archive>
    void load(Archive & ar, const unsigned int version)
    {
        ar & boost::serialization::base_object<base_t>(*this);
    }
    template <class Archive>
    void save(Archive & ar, const unsigned int version) const
    {
        ar & boost::serialization::base_object<base_t>(*this);
    }
    
    BOOST_SERIALIZATION_SPLIT_MEMBER()
};

template<class SymmGroup>
std::ostream& operator<<(std::ostream& os, DualIndex<SymmGroup> const & idx)
{
    os << "|";
    for (typename DualIndex<SymmGroup>::const_iterator it = idx.begin();
         it != idx.end();
         ++it)
    {
        os << "( " << (*it).lc << ","
                   << (*it).rc << ": "
                   << (*it).ls << "x"
                   << (*it).rs
           << " )";
    }
    os << "|";

    return os;
}

#endif
