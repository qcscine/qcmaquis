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

template<typename T, class SymmGroup>
class SymmetricMatrix
{
    SymmetricMatrix(DIndex<SymmGroup> rows_, DIndex<SymmGroup> cols_);
    
    DIndex<SymmGroup> rows() const;
    DIndex<SymmGroup> cols() const;
    
    // this will probably implement a similar interface to the standard Matrix
    // class...I would add methods as we go along. I don't know how we want to access
    // this thing, though...
};

// some example functions
template<typename T, class SymmGroup>
void gemm(
    SymmetricMatrix<T, SymmGroup> & A const,
    SymmetricMatrix<T, SymmGroup> & B const,
    SymmetricMatrix<T, SymmGroup> & C);

template<typename T, class SymmGroup>
void svd(
    SymmetricMatrix<T, SymmGroup> & M const,
    SymmetricMatrix<T, SymmGroup> &U, SymmetricMatrix<T, SymmGroup> &S, SymmetricMatrix<T, SymmGroup> &V,
    double truncation, DIndex<SymmGroup> maxdim,
    SVWhere where);

// a general reshape function
// The semantics of this are similar to the tensor class. In particular, I would
// assume that the names of indices in/out are the same, just re-ordered. The
// dimensions do not change, no rank changes -> no index fusion.

template<typename T, class SymmGroup, int R1, int R2, int R3>
void reshape(
    SymmetricMatrix<T, SymmGroup> & in const,
    boost::array<DIndex<SymmGroup>, R1> in_left,
    boost::array<DIndex<SymmGroup>, R2> in_right,
    SymmetricMatrix<T, SymmGroup> & out,
    boost::array<Index, R3> out_left,
    boost::array<Index, R1+R2-R3> out_right);
