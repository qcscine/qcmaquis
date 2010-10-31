#include <sstream>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/lambda/bind.hpp>

#include "indexing.h"
#include "symmetry.h"

namespace detail
{
    template<class A, class B> A first(std::pair<A, B> const & p) { return p.first; }
    template<class A, class B> A second(std::pair<A, B> const & p) { return p.second; }
    
    template<class SymmGroup> bool match(typename SymmGroup::charge c,
        std::pair<typename SymmGroup::charge, std::size_t> const & p)
    {
        return SymmGroup::SingletCharge == SymmGroup::fuse(c, p.first);
    }
}

template<class Matrix, class SymmGroup>
class block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
    
public:
    block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(), Index<SymmGroup> cols = Index<SymmGroup>())
    : rows_(rows), cols_(cols), data_(rows.size())
    , left_charge_map_(SymmGroup::get_map(rows.charges()))
    , right_charge_map_(SymmGroup::get_map(cols.charges()))
    {
        assert(rows_.size() == cols_.size());
        for (std::size_t k = 0; k < data_.size(); ++k)
            data_[k] = Matrix(rows_[k].second, cols_[k].second);
    }
    
    friend void swap(block_matrix<Matrix, SymmGroup> & x, block_matrix<Matrix, SymmGroup> & y)
    {
        std::swap(x.data_, y.data_);
        std::swap(x.rows_, y.rows_);
        std::swap(x.cols_, y.cols_);
        std::swap(x.left_charge_map_, y.left_charge_map_);
        std::swap(x.right_charge_map_, y.right_charge_map_);
    }
    
    block_matrix<Matrix, SymmGroup> & operator=(block_matrix<Matrix, SymmGroup> rhs)
    {
        swap(*this, rhs);
    }
    
    Index<SymmGroup> rows() const { return rows_; }
    Index<SymmGroup> cols() const { return cols_; }
    
    Matrix & atrow(charge c) { return data_[left_charge_map_[c]]; }
    Matrix const & atrow(charge c) const { return data_[left_charge_map_[c]]; }
    
    Matrix & atcol(charge c) { return data_[right_charge_map_[c]]; }
    Matrix const & atcol(charge c) const { return data_[right_charge_map_[c]]; }
    
    // use with caution!
    Matrix & operator[](std::size_t c) { return data_[c]; }
    Matrix const & operator[](std::size_t c) const { return data_[c]; }
    
    std::size_t n_blocks() const { return data_.size(); }
    
    friend void match_blocks(block_matrix<Matrix, SymmGroup> & a, block_matrix<Matrix, SymmGroup> & b)
    {
        // we match b to a
        for (std::size_t p1 = 0; p1 < a.n_blocks(); ++p1) {
            charge findme = a.cols()[p1].first;
            Index<SymmGroup> brows = b.rows(); // this way, we don't rely on rows() always returning a reference
            std::size_t p2 = std::find_if(brows.begin(), brows.end(),
                boost::lambda::bind(detail::match<SymmGroup>, findme, boost::lambda::_1)) - brows.begin();
            
            // now move p2 to p1 in the rhs
            if (p1 != p2) {
                iter_swap(b.data_.begin()+p2, b.data_.begin()+p1);
                iter_swap(b.rows_.begin()+p2, b.rows_.begin()+p1);
                iter_swap(b.cols_.begin()+p2, b.cols_.begin()+p1);
            }
            b.left_charge_map_ = SymmGroup::get_map(b.rows_.charges());
            b.right_charge_map_ = SymmGroup::get_map(b.cols_.charges());
        }
    }
    
    friend std::ostream& operator<<(std::ostream& os, block_matrix<Matrix, SymmGroup> const & m)
    {
        os << m.rows() << std::endl;
        os << m.cols() << std::endl;
        std::copy(m.data_.begin(), m.data_.end(), std::ostream_iterator<Matrix>(os, " "));
        os << std::endl;
        return os;
    }
    
    std::string description() const
    {
        std::ostringstream oss;
        oss << rows_ << cols_;
        return oss.str();
    }
    
protected:
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
    typename SymmGroup::map left_charge_map_, right_charge_map_;
};

// some example functions
template<class Matrix, class SymmGroup>
void gemm(
    block_matrix<Matrix, SymmGroup> & A,
    block_matrix<Matrix, SymmGroup> & B,
    block_matrix<Matrix, SymmGroup> & C)
{
    C = block_matrix<Matrix, SymmGroup>(A.rows(), B.cols());
    
    match_blocks(A, B);
    
    // Some checks
    // We should discuss whether we want to use asserts or something that doesn't disappear upon NDEBUG
    assert(A.n_blocks() == B.n_blocks() && B.n_blocks() == C.n_blocks());
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        // is the dimension of the sectors the same?
        assert(A.cols()[k].second == B.rows()[k].second);
        // is the charge of the column of the first the same as the row of the second?
        assert(SymmGroup::fuse(A.cols()[k].first, B.rows()[k].first) == SymmGroup::SingletCharge);
    }
    
    // this would require this function to be a friend of block_matrix<>
    // If we have gemm(A,B) -> C
    // std::transform(A.data_.begin(), A.data_.end(), B.data_.begin(), C.data_.begin(), gemm);
    // More likely:
    for (std::size_t k = 0; k < A.n_blocks(); ++k)
        gemm(A[k], B[k], C[k]);
}

// template<class Matrix, class SymmGroup>
// void svd(
//     block_matrix<Matrix, SymmGroup> & M,
//     block_matrix<Matrix, SymmGroup> &U, block_matrix<T, SymmGroup> &S, block_matrix<T, SymmGroup> &V,
//     double truncation, Index<SymmGroup> maxdim,
//     SVWhere where)
// {
//     /* basically analogous to gemm */
// }

// a general reshape function
// The semantics of this are similar to the tensor class. In particular, I would
// assume that the names of indices in/out are the same, just re-ordered. The
// dimensions do not change, no rank changes -> no index fusion.

template<class Matrix, class SymmGroup, int R1, int R2, int R3>
void reshape(
    block_matrix<Matrix, SymmGroup> & in,
    boost::array<Index<SymmGroup>, R1> in_left,
    boost::array<Index<SymmGroup>, R2> in_right,
    block_matrix<Matrix, SymmGroup> & out,
    boost::array<IndexName, R3> out_left,
    boost::array<IndexName, R1+R2-R3> out_right);
