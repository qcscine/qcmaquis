#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>

#include "indexing.h"
#include "symmetry.h"

namespace detail
{
    template<class SymmGroup, class Matrix>
    bool cmp0(
        boost::tuple<
            std::pair<typename SymmGroup::charge, std::size_t>,
            std::pair<typename SymmGroup::charge, std::size_t>,
            Matrix> const & a,
        boost::tuple<
            std::pair<typename SymmGroup::charge, std::size_t>,
            std::pair<typename SymmGroup::charge, std::size_t>,
            Matrix> const & b)
    {
        return boost::tuples::get<0>(a).first < boost::tuples::get<0>(b).first;
    }
    
    template<class SymmGroup, class Matrix>
    bool cmp1(
        boost::tuple<
            std::pair<typename SymmGroup::charge, std::size_t>,
            std::pair<typename SymmGroup::charge, std::size_t>,
            Matrix> const & a,
        boost::tuple<
            std::pair<typename SymmGroup::charge, std::size_t>,
            std::pair<typename SymmGroup::charge, std::size_t>,
            Matrix> const & b)
    {
        return boost::tuples::get<1>(a).first < boost::tuples::get<1>(b).first;
    }
}

template<class Matrix, class SymmGroup>
class block_matrix
{
public:
    block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(), Index<SymmGroup> cols = Index<SymmGroup>())
    : rows_(rows), cols_(cols), data_(rows.size())
    , left_charge_map_(SymmGroup::get_map(rows_.charges()))
    , right_charge_map_(SymmGroup::get_map(cols_.charges()))
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
    
    Matrix & atrow(typename SymmGroup::charge c) { return data_[left_charge_map_[c]]; }
    Matrix const & atrow(typename SymmGroup::charge c) const { return data_[left_charge_map_[c]]; }
    
    Matrix & atcol(typename SymmGroup::charge c) { return data_[right_charge_map_[c]]; }
    Matrix const & atcol(typename SymmGroup::charge c) const { return data_[right_charge_map_[c]]; }
    
    // use with caution!
    Matrix & operator[](std::size_t c) { return data_[c]; }
    Matrix const & operator[](std::size_t c) const { return data_[c]; }
    
    std::size_t n_blocks() const { return data_.size(); }
    
    void sort_rows()
    {
        // not readable, but elegant
        // FIXME: check whether this only uses swaps on the Matrix instead of explicit copying
        std::sort(boost::make_zip_iterator(rows_.begin(), cols_.begin(), data_.begin()),
            boost::make_zip_iterator(rows_.end(), cols_.end(), data_.end()),
            detail::cmp0<SymmGroup, Matrix>);
    }
    void sort_cols()
    {
         std::sort(boost::make_zip_iterator(rows_.begin(), cols_.begin(), data_.begin()),
                boost::make_zip_iterator(rows_.end(), cols_.end(), data_.end()),
                detail::cmp1<SymmGroup, Matrix>);
    }
    
protected:
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
    typename SymmGroup::charge_map left_charge_map_, right_charge_map_;
};

// some example functions
template<class Matrix, class SymmGroup>
void gemm(
    block_matrix<Matrix, SymmGroup> & A,
    block_matrix<Matrix, SymmGroup> & B,
    block_matrix<Matrix, SymmGroup> & C)
{
    C = block_matrix<Matrix, SymmGroup>(A.rows(), B.cols());
    
    A.sort_rows();
    B.sort_rows();
    C.sort_rows();
    
    // Some checks
    // We should discuss whether we want to use asserts or something that doesn't disappear upon NDEBUG
    assert(A.n_blocks() == B.n_blocks() && B.n_blocks() == C.n_blocks());
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        // is the dimension of the sectors the same?
        assert(A.cols().sectors[k].second == B.rows().sector[k].second);
        // is the charge of the column of the first the same as the row of the second?
        assert(A.cols().sectors[k].first == -B.rows().sector[k].first);
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
