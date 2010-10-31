// Symmetric Matrix, where symmetry of course means that the operator is the direct
// sum of irreducible representations of some internal symmetry group of the system.

// What's the difference between this and a tensor class? It's stupid!

// It doesn't know about its tensor structure. It only has a row and column
// dimension. It is therefore very easy to implement, at least without any
// fancy optimizations.
// It pushes a lot of bookkeeping to the user. In particular, all reshaping
// has to be carried out explicitly. The class does not have methods for this,
// but there are global functions. It is probably easiest to have a general one,
// which I can implement basically by copying to/from my tensor class, and then
// special ones for the reshapes usually needed, in particular the
// (alpha,sigma),beta -> alpha,(beta,sigma) reshape
// We can write these one after another when we see what shows up as most
// costly operations in benchmarks. This can be done in a 'non-intrusive' way,
// i.e. no hacking around template specializations, hooks, or other intransparent stuff.

// The huge advantage of this approach is that it makes the structure of the underlying
// operations very transparent -> easier to optimize, easier to factorize into
// smaller parts.

template<class Matrix, class SymmGroup>
class block_matrix
{
public:
    block_matrix(DIndex<SymmGroup> rows = DIndex<SymmGroup>(), DIndex<SymmGroup> cols = DIndex<SymmGroup>())
    {
        /* to be implemented */
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
    
    DIndex<SymmGroup> rows() const { return rows_; }
    DIndex<SymmGroup> cols() const { return cols_; }
    
    // Since these matrices are always block-diagonal, at least in the sense of uniquely
    // mapping rows charge to column charge, only one index (row/col) is necessary to index
    // them. I therefore would not provide an operator(), since it is not immediately clear
    // how that would behave.
    
    Matrix & atrow(typename SymmGroup::charge c) { return data_[left_charge_map[c]]; }
    Matrix const & atrow(typename SymmGroup::charge, typename SymmGroup::charge) const { return data_[left_charge_map[c]]; }
    
    Matrix & atcol(typename SymmGroup::charge) { return data_[right_charge_map[c]]; }
    Matrix const & atcol(typename SymmGroup::charge, typename SymmGroup::charge) const { return data_[right_charge_map[c]]; }
    
    // use with caution!
    Matrix & operator[](std::size_t c) { return data_[c]; }
    Matrix const & operator[](std::size_t c) const { return data_[c]; }
    
    std::size_t n_blocks() const { return data_.size(); }
    
    void sort_rows();
    void sort_cols();
    
protected:
    std::vector<Matrix> data_;
    DIndex<SymmGroup> rows_, cols_;
    typename SymmGroup::charge_map left_charge_map_, right_charge_map_;
};

// some example functions
template<class Matrix, class SymmGroup>
void gemm(
    block_matrix<Matrix, SymmGroup> & A const,
    block_matrix<Matrix, SymmGroup> & B const,
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

template<class Matrix, class SymmGroup>
void svd(
    block_matrix<Matrix, SymmGroup> & M const,
    block_matrix<Matrix, SymmGroup> &U, block_matrix<T, SymmGroup> &S, block_matrix<T, SymmGroup> &V,
    double truncation, DIndex<SymmGroup> maxdim,
    SVWhere where)
{
    /* basically analogous to gemm */
}

// a general reshape function
// The semantics of this are similar to the tensor class. In particular, I would
// assume that the names of indices in/out are the same, just re-ordered. The
// dimensions do not change, no rank changes -> no index fusion.

template<class Matrix, class SymmGroup, int R1, int R2, int R3>
void reshape(
    block_matrix<Matrix, SymmGroup> & in const,
    boost::array<DIndex<SymmGroup>, R1> in_left,
    boost::array<DIndex<SymmGroup>, R2> in_right,
    block_matrix<Matrix, SymmGroup> & out,
    boost::array<Index, R3> out_left,
    boost::array<Index, R1+R2-R3> out_right);
