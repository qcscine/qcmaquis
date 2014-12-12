/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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
    class QnBlock
    {
        typedef typename SymmGroup::charge charge;

    public:
        QnBlock() {}
        QnBlock(charge lc_, charge rc_, std::size_t ls_, std::size_t rs_)
            : lc(lc_), rc(rc_), ls(ls_), rs(rs_) {}

        bool operator==(QnBlock const & o) const
        {
            return lc == o.lc && rc == o.rc && ls == o.ls && rs == o.rs;
        }

        typename SymmGroup::charge lc;
        typename SymmGroup::charge rc;
        std::size_t                ls;
        std::size_t                rs;
    };

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
{
    typedef std::vector<dual_index_detail::QnBlock<SymmGroup> > data_type;
    
public:
    typedef typename SymmGroup::charge charge;
    typedef typename data_type::value_type value_type;
    
    typedef typename data_type::iterator iterator;
    typedef typename data_type::const_iterator const_iterator;
    
    typedef typename data_type::reverse_iterator reverse_iterator;
    typedef typename data_type::const_reverse_iterator const_reverse_iterator;
    
    typedef basis_iterator_<SymmGroup> basis_iterator;
    
    DualIndex() : sorted_(true) {}
    
    std::size_t left_block_size(charge r, charge c) const {
        std::size_t pos = position(dual_index_detail::QnBlock<SymmGroup>(r,c,0,0));
        return (*this)[pos].ls;
    }
    std::size_t right_block_size(charge r, charge c) const {
        std::size_t pos = position(dual_index_detail::QnBlock<SymmGroup>(r,c,0,0));
        return (*this)[pos].rs;
    }
    
    std::size_t position(charge row, charge col) const
    {
        const_iterator match;
        if (sorted_)
            match = std::lower_bound(data_.begin(), data_.end(), dual_index_detail::QnBlock<SymmGroup>(row,col,0,0), dual_index_detail::gt<SymmGroup>());
        else
            match = std::find_if(data_.begin(), data_.end(), dual_index_detail::is_first_equal<SymmGroup>(row,col));
        
        if (match != data_.end() && ((*match).lc != row || (*match).rc != col))
            match = data_.end();
        return std::distance(data_.begin(), match);
    }

    bool has(charge row, charge col) const
    {
        if (sorted_)
            return std::binary_search(data_.begin(), data_.end(), dual_index_detail::QnBlock<SymmGroup>(row,col,0,0), dual_index_detail::gt<SymmGroup>());
        else
            return std::find_if(data_.begin(), data_.end(),
                                dual_index_detail::is_first_equal<SymmGroup>(row,col)) != data_.end();
    }

    bool is_sorted() const { return sorted_; }
    
    void sort()
    {
        std::sort(data_.begin(), data_.end(), dual_index_detail::gt<SymmGroup>());
        sorted_ = true;
    }
    
    std::size_t insert(value_type const & x)
    {
        if (sorted_) {
            std::size_t d = destination(x);
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
        {
            (*this)[k].lc = SymmGroup::fuse((*this)[k].lc, diff);
            (*this)[k].rc = SymmGroup::fuse((*this)[k].rc, diff);
        }
    }
    
    bool operator==(DualIndex const & o) const
    {
        return (data_.size() == o.size()) && std::equal(data_.begin(), data_.end(), o.begin());
    }

    bool operator!=(DualIndex const & o) const
    {
        return !( *this == o );
    }

    basis_iterator basis_begin() const
    {
        assert( data_.size() > 0 );
        return basis_iterator(*this);
    }
    
    std::size_t sum_of_left_sizes() const
    {
        return std::accumulate(data_.begin(), data_.end(), 0,
                               boost::lambda::_1
                               + boost::lambda::bind(&dual_index_detail::QnBlock<SymmGroup>::ls, boost::lambda::_2)
                              );
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

    charge & left_charge(std::size_t k) { return data_[k].lc; }
    charge & right_charge(std::size_t k) { return data_[k].rc; }
    std::size_t & left_size(std::size_t k) { return data_[k].ls; }
    std::size_t & right_size(std::size_t k) { return data_[k].rs; }

    charge const & left_charge(std::size_t k) const { return data_[k].lc; }
    charge const & right_charge(std::size_t k) const { return data_[k].rc; }
    std::size_t const & left_size(std::size_t k) const { return data_[k].ls; }
    std::size_t const & right_size(std::size_t k) const { return data_[k].rs; }

    void resize(std::size_t sz) { data_.resize(sz); }
    
    value_type & operator[](std::size_t p) { return data_[p]; }
    value_type const & operator[](std::size_t p) const { return data_[p]; }
    
    std::size_t size() const { return data_.size(); }
    
    iterator erase(iterator p) { iterator r = data_.erase(p); return r; }
    iterator erase(iterator a, iterator b) { iterator r = data_.erase(a,b); return r; }

    friend void swap(DualIndex & a, DualIndex & b)
    {
        using std::swap;
        swap(a.data_,   b.data_);
        swap(a.sorted_, b.sorted_);
    }
    
private:
    data_type data_;
    bool sorted_;
    
    void push_back(value_type const & x){
        data_.push_back(x);
    }
    
    std::size_t destination(value_type const & x) const
    {
        return std::find_if(data_.begin(), data_.end(),
                            boost::lambda::bind(dual_index_detail::lt<SymmGroup>,
                                                boost::lambda::_1,
                                                x)) - data_.begin();
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
        ar["DualIndex"] >> data_;
    }
    template <class Archive>
    void save(Archive & ar) const
    {
        ar["DualIndex"] << data_;
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
