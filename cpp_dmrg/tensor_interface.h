#ifndef TENSOR_INTERFACE_H
#define TENSOR_INTERFACE_H

namespace tensor {
    // wrap a tensor into a matrix for all sorts of lazy evaluation functions
    // replaces all the transposes/reshapes in Numpy
    // Conjugation or other element-wise operations could be stored in this
    // The template parameter allows inlining. Do we believe that this is necessary?
    struct NoOp { template<class T> static T operator()(T v) { return v; } };
    struct ConjOp { template<class T> static T operator()(T, v) { return conj(v); } };
    
    template<class Op> struct conj_trait { };
    template<> struct conj_trait<NoOp> { typedef ConjOp result; };
    template<> struct conj_trait<ConfOp> { typedef NoOp result; };

    // ElementOp is purely a flag, a converting constructor will be required
    template<typename T, unsigned int R1, unsigned int R2, class SymmGroup, class ElementOp = NoOp>
    class tensor_matrix_wrapper;
    
    template<typename T, unsigned int R1, unsigned int R2, class SymmGroup, class ElementOpIn>
    tensor_matrix_wrapper<T, R1, R2, SymmGroup, typename conj_trait<ElementOpIn>::result>
    conj(tensor_matrix_wrapper<T, R1, R2, SymmGroup, ElementOpIn> const & v) { return v; }

    template<typename T, unsigned int R, class SymmGroup>
    class tensor
    {
    public:
        tensor();
        tensor(const tensor<T, R, SymmGroup>&);
        tensor(boost::array<Index, R> names, boost::array<DIndex, R> dimensions);
    
        DIndex<SymmGroup> get_dimensions(Index);
        
        template<int RR> void rename(boost::array<Index, RR> from, boost::array<Index, RR> to);
        void rename(Index from, Index to);
    
        // Input for matrix-based operations
        template<int R1, int R2>
        tensor_matrix_wrapper<T, R1, R2, SymmGroup> operator()(boost::array<Index, R1>, boost::array<Index, R2>);
    
        // Random access: cumbersome because of two indices
        // Any better ideas?
        T & operator()(boost::array<typename SymmGroup::charge, R> const &, boost::array<std::size_t, R> const &)
    };

    // Multiply the first two tensors into the third one
    // Internally, perform conversion to matrix form and use callbacks to optimized functions
    // Why is no parameter const? Because of potential reshaping, which would require us to make basically
    // all members mutable
    template<int R1, int R2, int R3, typename T, class SymmGroup>
    void mult(
        tensor_matrix_wrapper<T, R1, R2, SymmGroup> &,
        tensor_matrix_wrapper<T, R2, R3, SymmGroup> &,
        tensor<T, R1+R3, SymmGroup> &);
    
    /* SVD tensor into three tensors, one of which is diagonal
       Params:
       @M: to be decomposed
       @U,V: take a wild guess
       @S: why a pointer? maybe we don't care about singular values, in that case pass NULL
       @name: the 'internal' indices on U,V, need a name
       @truncation: where to truncate singular values
       @where: possible return U*S, S*V, or U*sqrt(S) sqrt(S)*V
    */
    
    /* SVD Design question: the below function doesn't allow to estimate what was discarded
     in the truncation. This may be computationally advantageous, since only the actually
     needed singular values/vectors have to be calculated. We don't have an implementation
     that exploits this, though. */
    
    enum SVWhere { SV_TO_LEFT, SV_TO_RIGHT, SV_TO_NONE };
    
    template<int R1, int R2, typename T, class SymmGroup>
    void svd(
        tensor_matrix_wrapper<T, R1, R2, SymmGroup> &M,
        tensor<T, R1+1, SymmGroup> &U,
        tensor<T, R2+2, SymmGroup> &V,
        tensor<T, 2, SymmGroup> *S,
        Index name,
        double truncation,
        DIndex<SymmGroup> maxdim,
        SVWhere where);
    
    template<int R1, int R2, typename T, class SymmGroup>
    void qr(
        tensor_matrix_wrapper<T, R1, R2, SymmGroup> &M,
        tensor<T, R1+1, SymmGroup> &Q,
        tensor<T, R2+1, SymmGroup> &R);
    
} // namespace tensor

#endif
