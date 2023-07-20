/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TENSOR_DUAL_INDEX_H
#define TENSOR_DUAL_INDEX_H

#include "dmrg/block_matrix/indexing_stable.hpp"
#include <vector>

#include <boost/tuple/tuple_comparison.hpp>
#include <boost/tuple/tuple.hpp>
#include <boost/serialization/nvp.hpp>

namespace dual_index_detail
{
    template <class SymmGroup>
    class QnBlock
    {
        using charge = typename SymmGroup::charge;

    public:
        QnBlock() = default;
        QnBlock(charge lc_, charge rc_, std::size_t ls_, std::size_t rs_)
            : lc(lc_), rc(rc_), ls(ls_), rs(rs_) {}

        bool operator==(QnBlock const & o) const
        {
            return lc == o.lc && rc == o.rc && ls == o.ls && rs == o.rs;
        }

        charge lc;
        charge rc;
        std::size_t                ls{};
        std::size_t                rs{};
    };

    template<class SymmGroup>
    struct gt {
        bool operator()(QnBlock<SymmGroup> const & a,
                        QnBlock<SymmGroup> const & b)
        {
            if (a.lc > b.lc) {
                return true;
            } else if (a.lc < b.lc) {
                return false;
            } else {
                return a.rc > b.rc;
            }
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
        if (a.lc < b.lc) {
            return true;
        } else if (a.lc > b.lc) {
            return false;
        } else {
            return a.rc < b.rc;
        }
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

    template<class SymmGroup>
    class is_first_equal_row
    {
    public:
        is_first_equal_row(typename SymmGroup::charge c1) : c1_(c1) { }

        bool operator()(QnBlock<SymmGroup> const & x) const
        {
            return x.lc == c1_;
        }

    private:
        typename SymmGroup::charge c1_;
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
    using data_type = std::vector<dual_index_detail::QnBlock<SymmGroup>>;
    
public:
    using charge = typename SymmGroup::charge;
    using value_type = typename data_type::value_type;
    
    using iterator = typename data_type::iterator;
    using const_iterator = typename data_type::const_iterator;
    
    using reverse_iterator = typename data_type::reverse_iterator;
    using const_reverse_iterator = typename data_type::const_reverse_iterator;
    
    using basis_iterator = basis_iterator_<SymmGroup>;
    
    DualIndex() : sorted_(true) {}
    
    std::size_t left_block_size(charge r, charge c) const {
        std::size_t pos = position(value_type(r,c,0,0));
        return (*this)[pos].ls;
    }
    std::size_t right_block_size(charge r, charge c) const {
        std::size_t pos = position(value_type(r,c,0,0));
        return (*this)[pos].rs;
    }
    
    std::size_t position(charge row, charge col) const
    {
        const_iterator match;
        if (sorted_) {
            match = std::lower_bound(data_.begin(), data_.end(), value_type(row,col,0,0), dual_index_detail::gt<SymmGroup>());
        } else {
            match = std::find_if(data_.begin(), data_.end(), dual_index_detail::is_first_equal<SymmGroup>(row,col));
        }
        
        if (match != data_.end() && ((*match).lc != row || (*match).rc != col)) {
            match = data_.end();
        }
        return std::distance(data_.begin(), match);
    }

    bool has(charge row, charge col) const
    {
        if (sorted_) {
            return std::binary_search(data_.begin(), data_.end(), value_type(row,col,0,0), dual_index_detail::gt<SymmGroup>());
        } else {
            return std::find_if(data_.begin(), data_.end(),
                                dual_index_detail::is_first_equal<SymmGroup>(row,col)) != data_.end();
        }
    }

    const_iterator left_lower_bound(charge row) const
    {
        if (sorted_) {
            return std::lower_bound(
                data_.begin(), data_.end(),
                value_type(row, SymmGroup::IdentityCharge,0,0),
                dual_index_detail::gt_row<SymmGroup>());
        }
        else {
            return std::find_if(
                data_.begin(), data_.end(),
                dual_index_detail::is_first_equal_row<SymmGroup>(row));
        }
    }

    std::pair<const_iterator, const_iterator> left_equal_range(charge row) const
    {
        if (sorted_) {
            return std::equal_range(
                data_.begin(), data_.end(),
                value_type(row, SymmGroup::IdentityCharge,0,0),
                dual_index_detail::gt_row<SymmGroup>());
        }
        else {
          throw std::runtime_error("Not implemented for unsorted");
//            return std::make_pair(std::find_if(
 //               data_.begin(), data_.end(),
  //              dual_index_detail::is_first_equal_row<SymmGroup>(row)),
   //                std::find_if(
    //            data_.begin(), data_.end(),
     //           !(dual_index_detail::is_first_equal_row<SymmGroup>(row))));
        }
    }

    bool left_has(charge row) const
    {
        const_iterator it = left_lower_bound(row);
        return (it != data_.end() && row == it->lc);
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
            [&](const auto& acc, const auto& x){
                return acc + x.ls;
            });
    }

    std::size_t memory_size() const
    {
        return std::accumulate(data_.begin(), data_.end(), 0,
            [](const auto& acc, const auto& x){
                return acc + x.ls * x.rs;
            });
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
      return std::distance(
          data_.begin(),
            std::find_if(
              data_.begin(), data_.end(), 
              [x](const auto& e){
                  return dual_index_detail::lt<SymmGroup>(e, x);
              }));
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
