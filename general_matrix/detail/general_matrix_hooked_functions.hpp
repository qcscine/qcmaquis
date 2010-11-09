#include <boost/numeric/bindings/lapack/driver/gesdd.hpp>

namespace blas {
    template <typename T, typename MemoryBlock>
        class general_matrix;
}

//
// Hooked general matrix functions
//
namespace blas { namespace detail {

#define MATRIX_MATRIX_MULTIPLY(T) \
    template <typename MemoryBlock> \
    const general_matrix<T,MemoryBlock> matrix_matrix_multiply(general_matrix<T,MemoryBlock> const& lhs, general_matrix<T,MemoryBlock> const& rhs) \
    { \
        assert( lhs.num_columns() == rhs.num_rows() ); \
        general_matrix<T,MemoryBlock> result(lhs.num_rows(),rhs.num_columns()); \
        boost::numeric::bindings::blas::gemm \
            ( \
               typename general_matrix<T,MemoryBlock>::value_type(1), \
               lhs, \
               rhs, \
               typename general_matrix<T,MemoryBlock>::value_type(1), \
               result \
            ); \
        return result; \
    }
IMPLEMENT_FOR_ALL_BLAS_TYPES(MATRIX_MATRIX_MULTIPLY)
#undef MATRIX_MATRIX_MULTIPLY

    template <typename T, typename MemoryBlock>
    const general_matrix<T,MemoryBlock> matrix_matrix_multiply(general_matrix<T,MemoryBlock> const& lhs, general_matrix<T,MemoryBlock> const& rhs)
    {
        assert( lhs.num_columns() == rhs.num_rows() );

        // Simple matrix matrix multiplication
        general_matrix<T,MemoryBlock> result(lhs.num_rows(),rhs.num_columns());
        for(std::size_t i=0; i < lhs.num_rows(); ++i)
        {
            for(std::size_t j=0; j<rhs.num_columns(); ++j)
            {
                for(std::size_t k=0; k<lhs.num_columns(); ++k)
                {
                        result(i,j) += lhs(i,k) * rhs(k,j);
                }
            }
        }
        return result;
    } 

// This seems to be the best solution for the *_ASSIGN dispatchers at the moment even though they call functions within the detail namespace
#define PLUS_MINUS_ASSIGN(T) \
    template <typename MemoryBlock> \
    void plus_and_minus_assign_impl(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs, typename general_matrix<T,MemoryBlock>::value_type const& sign) \
    { \
        assert( m.num_columns() == rhs.num_columns() && m.num_rows() == rhs.num_rows() ); \
        if(!(m.is_shrinkable() || rhs.is_shrinkable()) ) \
        { \
            boost::numeric::bindings::blas::detail::axpy( m.num_rows() * m.num_columns(), sign, &(*rhs.column(0).first), 1, &(*m.column(0).first), 1); \
        } \
        else \
        { \
            for(std::size_t j=0; j < m.num_columns(); ++j) \
                boost::numeric::bindings::blas::detail::axpy( m.num_rows(), sign, &(*rhs.column(j).first), 1, &(*m.column(j).first), 1); \
        } \
    } \
    template <typename MemoryBlock> \
    void plus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs) \
        { plus_and_minus_assign_impl(m, rhs, typename general_matrix<T,MemoryBlock>::value_type(1)); } \
    template <typename MemoryBlock> \
    void minus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs) \
        { plus_and_minus_assign_impl(m, rhs, typename general_matrix<T,MemoryBlock>::value_type(-1)); }
IMPLEMENT_FOR_ALL_BLAS_TYPES(PLUS_MINUS_ASSIGN)
#undef PLUS_MINUS_ASSIGN

    template <typename T,typename MemoryBlock>
    void plus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs)
    {
        m.plus_assign(rhs);
    }

    template <typename T, typename MemoryBlock>
    void minus_assign(general_matrix<T,MemoryBlock>& m, general_matrix<T,MemoryBlock> const& rhs)
    {
        m.minus_assign(rhs);
    }


#define MULTIPLIES_ASSIGN(T) \
    template <typename MemoryBlock> \
    void multiplies_assign(general_matrix<T,MemoryBlock>& m, T const& t) \
    { \
        if( !(m.is_shrinkable()) ) \
        { \
            boost::numeric::bindings::blas::detail::scal( m.num_rows()*m.num_columns(), t, &(*m.column(0).first), 1 ); \
        } \
        else \
        { \
            for(std::size_t j=0; j <m.num_columns(); ++j) \
                boost::numeric::bindings::blas::detail::scal( m.num_rows(), t, &(*m.column(j).first), 1 ); \
        } \
    }
    IMPLEMENT_FOR_ALL_BLAS_TYPES(MULTIPLIES_ASSIGN)
#undef MULTIPLIES_ASSIGN


    template <typename T, typename MemoryBlock>
    void multiplies_assign(general_matrix<T,MemoryBlock>& m, T const& t)
    {
        m.multiplies_assign(t);
    }

}} // namespace detail, namespace blas

