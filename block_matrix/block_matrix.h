#include <sstream>

#include <boost/tuple/tuple.hpp>
#include <boost/iterator/zip_iterator.hpp>
#include <boost/lambda/bind.hpp>

#include "indexing.h"
#include "symmetry.h"

template<class Matrix, class SymmGroup>
class block_matrix
{
private:
    typedef typename SymmGroup::charge charge;
    
public:
    block_matrix(Index<SymmGroup> rows = Index<SymmGroup>(), Index<SymmGroup> cols = Index<SymmGroup>())
    : rows_(rows)
    , cols_(cols)
    {
        assert(rows_.size() == cols_.size());
        for (std::size_t k = 0; k < rows_.size(); ++k)
            data_.push_back(Matrix(rows_[k].second, cols_[k].second));
    }
    
    block_matrix(charge c, Matrix const & m)
    {
        rows_.push_back(std::make_pair(c, m.num_rows()));
        cols_.push_back(std::make_pair(c, m.num_cols()));
        data_.push_back(m);
    }
    
    block_matrix & operator+=(block_matrix const & rhs)
    {
        for (std::size_t k = 0; k < rhs.n_blocks(); ++k)
        {
            charge rhs_rc = rhs.rows_[k].first;
            std::size_t goesto = rows_.at(rhs_rc);
            if (goesto == rows_.size()) { // it's a new block
                std::size_t i1 = rows_.insert(rhs.rows_[k]);
                std::size_t i2 = cols_.insert(rhs.cols_[k]);
                assert(i1 == i2);
                data_.insert(data_.begin() + i1, rhs.data_[k]);
            } else { // this block exists already -> pass to Matrix
                data_[goesto] += rhs.data_[k];
            }
        }
        return *this;
    }
    
    void insert_block(boost::tuple<Matrix, charge, charge> const & block)
    {
        Matrix const & mtx = boost::tuples::get<0>(block);
        
        std::pair<charge, std::size_t>
        p1 = std::make_pair(boost::tuples::get<1>(block), mtx.num_rows()),
        p2 = std::make_pair(boost::tuples::get<2>(block), mtx.num_columns());
        
        std::size_t i1 = rows_.insert(p1);
        std::size_t i2 = cols_.insert(p2);
        assert(i1 == i2);
        data_.insert(data_.begin() + i1, mtx);
    }
    
    friend void swap(block_matrix & x, block_matrix & y)
    {
        std::swap(x.data_, y.data_);
        std::swap(x.rows_, y.rows_);
        std::swap(x.cols_, y.cols_);
    }
    
    block_matrix & operator=(block_matrix rhs)
    {
        swap(*this, rhs);
        return *this;
    }
    
    Index<SymmGroup> left_basis() const { return rows_; }
    Index<SymmGroup> right_basis() const { return cols_; }
    
    std::size_t n_blocks() const { return data_.size(); }
    
    std::string description() const
    {
        std::ostringstream oss;
        oss << rows_ << cols_;
        return oss.str();
    }
    
    Matrix & operator[](std::size_t c) { return data_[c]; }
    Matrix const & operator[](std::size_t c) const { return data_[c]; }
private:
    
    std::vector<Matrix> data_;
    Index<SymmGroup> rows_, cols_;
};    

template<class Matrix, class SymmGroup>    
std::ostream& operator<<(std::ostream& os, block_matrix<Matrix, SymmGroup> const & m)
{
    os << "Left HS: " << std::endl << m.left_basis();
    os << "Right HS: " << std::endl << m.right_basis();
//    for (std::size_t k = 0; k < m.n_blocks(); ++k)
//        os << "Block (" << m.num_rows()[k].first << "," << m.num_cols()[k].first << "):" << std::endl << m[k];
    os << std::endl;
    return os;
}

// some example functions
template<class Matrix, class SymmGroup>
void gemm(block_matrix<Matrix, SymmGroup> & A,
          block_matrix<Matrix, SymmGroup> & B,
          block_matrix<Matrix, SymmGroup> & C)
{
    C = block_matrix<Matrix, SymmGroup>(A.left_basis(), B.right_basis());
    
    // Some checks
    // We should discuss whether we want to use asserts or something that doesn't disappear upon NDEBUG
    assert(A.n_blocks() == B.n_blocks() && B.n_blocks() == C.n_blocks());
    for (std::size_t k = 0; k < A.n_blocks(); ++k) {
        // is the charge of the column of the first the same as the row of the second?
        assert(A.right_basis()[k].first == B.left_basis()[k].first);
        // is the dimension of the sectors the same?
        assert(A.right_basis()[k].second == B.left_basis()[k].second);
    }
    
    for (std::size_t k = 0; k < A.n_blocks(); ++k)
        C[k] = A[k] * B[k];
}

template<class Matrix, class SymmGroup>
void svd(block_matrix<Matrix, SymmGroup> const & M,
         block_matrix<Matrix, SymmGroup> & U,
         block_matrix<Matrix, SymmGroup> & V,
         block_matrix<Matrix, SymmGroup> & S)
{
    Index<SymmGroup> r = M.left_basis(), c = M.right_basis(), m = M.left_basis();
    for (std::size_t i = 0; i < M.n_blocks(); ++i)
        m[i].second = std::min(r[i].second, c[i].second);
    
    U = block_matrix<Matrix, SymmGroup>(r, transpose(m));
    V = block_matrix<Matrix, SymmGroup>(m, c);
    S = block_matrix<Matrix, SymmGroup>(m, transpose(m));
    
    for (std::size_t k = 0; k < M.n_blocks(); ++k)
    {
        std::vector<double> sk;
        svd(M[k], U[k], V[k], sk);
        for (std::size_t l = 0; l < sk.size(); ++l)
            S[k](l,l) = sk[l];
    }
}

template<class Matrix, class SymmGroup>
block_matrix<Matrix, SymmGroup> reshape_right(block_matrix<Matrix, SymmGroup>, Index<SymmGroup> phys_space);
