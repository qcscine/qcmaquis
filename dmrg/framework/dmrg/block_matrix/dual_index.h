/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef TENSOR_DUAL_INDEX_H
#define TENSOR_DUAL_INDEX_H

#include "dmrg/block_matrix/indexing_stable.hpp"
#include <vector>
#include <functional>

#include <boost/serialization/nvp.hpp>

namespace dual_index_detail
{
    /**
     * @brief QnBlock contains the information for a quantum number block.
     *        This includes left, right charges and left, right sizes.
     */
    template <class SymmGroup>
    class QnBlock {
        using charge = typename SymmGroup::charge;

    public:
        QnBlock() = default;
        QnBlock(charge lc_, charge rc_, std::size_t ls_, std::size_t rs_)
            : lc(lc_), rc(rc_), ls(ls_), rs(rs_) {}

        bool operator==(const QnBlock& other) const {
            return lc == other.lc && rc == other.rc && ls == other.ls && rs == other.rs;
        }
        // Left charge takes precedance over right charge
        bool operator>(const QnBlock& other) const {
          if (lc == other.lc) {
            return rc > other.rc;
          } else {
            return lc > other.lc;
          }
        }
        bool operator<(const QnBlock& other) const {
          return other > *this;
        }
        bool charges_equal(const QnBlock& other) const {
          return lc == other.lc && rc == other.rc;
        }

        charge lc;          // left charge
        charge rc;          // right charge
        std::size_t ls{};   // left size
        std::size_t rs{};   // right size
    };

    // --- Convenience functions ---
    const auto is_left_charge_greater = [](const auto& QnBlockA, const auto& QnBlockB) -> bool {
      return QnBlockA.lc > QnBlockB.lc;
    };

    const auto are_charges_equal = [](const auto& QnBlock, const auto& left_charge, const auto& right_charge) -> bool {
      return QnBlock.lc == left_charge && QnBlock.rc == right_charge;
    };

    /** @brief This is used in case the index is not sorted.
     *         Checks only for equality of the left charge.
     */
    const auto is_left_charge_equal = [](const auto& QnBlock, const auto& charge) -> bool {
      return QnBlock.lc == charge;
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


/**
 * @brief DualIndex stores the information about the quantum number blocks
 *        of a block matrix. This is achieved by storing a vector of QnBlock
 *        where each element corresponds to a specific block of the block matrix.
 * @warning  The DualIndex is sorted in descending order.
 */
template<class SymmGroup>
class DualIndex
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

    /**
     * @brief Returns the postion of the element with charge equal to (row, col)
     *        if exists, otherwise returns the position of the last element.
     */
    std::size_t position(charge row, charge col) const
    {
        const_iterator match;
        if (sorted_) {
            // Finds first element whos charge is not greater than (row, col)
            // Note that the DualIndex is sorted in descending order.
            match = std::lower_bound(
                data_.begin(),
                data_.end(),
                value_type(row,col,0,0),
                std::greater{});
        } else {
            match = std::find_if(
                data_.begin(),
                data_.end(),
                [&row, &col](const auto& QnBlock){
                   return dual_index_detail::are_charges_equal(QnBlock, row, col);
                });
        }
        
        // If the element is not found, return the position of the last element.
        bool not_found = match != data_.end() && ((*match).lc != row || (*match).rc != col);
        if (not_found) {
            match = data_.end();
        }
        return std::distance(data_.begin(), match);
    }

    /** @brief Checks if quantum number matching (row, col) exists. */
    bool has(charge row, charge col) const
    {
        if (sorted_) {
            return std::binary_search(
                data_.begin(),
                data_.end(),
                value_type(row,col,0,0),
                std::greater{});
        } else {
            return std::find_if(
                data_.begin(),
                data_.end(),
                [&row, &col](const auto& QnBlock){
                   return dual_index_detail::are_charges_equal(QnBlock, row, col);
                }) != data_.end();
        }
    }

    /** @brief Finds first element matching with left charge equal to row. */
    const_iterator left_lower_bound(charge row) const
    {
        if (sorted_) {
            return std::lower_bound(
                data_.begin(), data_.end(),
                value_type(row, SymmGroup::IdentityCharge,0,0),
                dual_index_detail::is_left_charge_greater);
        }
        else {
            return std::find_if(
                data_.begin(), data_.end(),
                [&row](const auto& QnBlock){
                   return dual_index_detail::is_left_charge_equal(QnBlock, row);
                });
        }
    }

    /** @brief Returns all the elements with the same left charge as row. */
    std::pair<const_iterator, const_iterator> left_equal_range(charge row) const
    {
        if (sorted_) {
            return std::equal_range(
                data_.begin(), data_.end(),
                value_type(row, SymmGroup::IdentityCharge,0,0),
                dual_index_detail::is_left_charge_greater);
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

    /** @brief Checks if block with charge equal to row exists. */
    bool left_has(charge row) const
    {
        const_iterator it = left_lower_bound(row);
        return (it != data_.end() && row == it->lc);
    }

    /** @brief Checks if sorted flag is set */
    bool is_sorted() const { return sorted_; }
    
    /** @brief Sorts the DualIndex in descending order. */
    void sort()
    {
        // kszenes: Why are we sorting in descending order?
        std::sort(data_.begin(), data_.end(), std::greater{});
        sorted_ = true;
    }
    
    /** @brief Inserts an element into the DualIndex and returns its position. */
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
    
    /** @brief Shifts bothe left and right charge by diff. */
    void shift(charge diff)
    {
        for (auto& x : data_) {
          x.lc = SymmGroup::fuse(x.lc, diff);
          x.rc = SymmGroup::fuse(x.rc, diff);
        }
    }
    
    bool operator==(DualIndex const & o) const
    {
        return (data_.size() == o.size())
          && std::equal(data_.begin(), data_.end(), o.begin());
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
    
    /** @brief Computes sum or rows in block_matrix. */
    std::size_t sum_of_left_sizes() const
    {
        return std::accumulate(data_.begin(), data_.end(), 0,
            [&](const auto& acc, const auto& x){
                return acc + x.ls;
            });
    }

    /** @brief Computes total number of elements in block_matrix. */
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
    charge const & left_charge(std::size_t k) const { return data_[k].lc; }
    charge const & right_charge(std::size_t k) const { return data_[k].rc; }
    std::pair<charge, charge> & charges(std::size_t k) const {
      return {left_charge(k), right_charge(k)};
    }
    std::size_t & left_size(std::size_t k) { return data_[k].ls; }
    std::size_t & right_size(std::size_t k) { return data_[k].rs; }
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
          std::find_if(data_.begin(),data_.end(),
            [x](const auto& e){ return e < x; }));
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
std::ostream& operator<<(std::ostream& os, const dual_index_detail::QnBlock<SymmGroup>& block)
{
      os << "( " << block.lc << ","
                 << block.rc << ": "
                 << block.ls << "x"
                 << block.rs
         << " )";
    return os;
}

template<class SymmGroup>
std::ostream& operator<<(std::ostream& os, const DualIndex<SymmGroup>& idx)
{
    os << "|";
    for (const auto& QnBlock : idx) {
      os << QnBlock;

    }
    os << "|";

    return os;
}

#endif
