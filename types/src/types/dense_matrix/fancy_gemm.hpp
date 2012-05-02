/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef FANCY_GEMM_HPP
#define FANCY_GEMM_HPP

#include <boost/numeric/bindings/tag.hpp>
#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>

#include "types/utils/traits.hpp"

namespace blas {
    template <typename T, typename MemoryBlock>
    class dense_matrix;
}

namespace blas {
    
    namespace detail {
        template<class T, class MemoryBlock, class Tag>
        struct TransConjWrapper
        {
            TransConjWrapper(dense_matrix<T, MemoryBlock> const & ref) : ptr(&ref) { }
            dense_matrix<T, MemoryBlock> const * ptr;
        };
        
        template<class Tag>
        struct resolve_to_ublas;
        
        template<> struct resolve_to_ublas<NoTranspose>
        { typedef boost::numeric::bindings::tag::column_major value; };
        
        template<> struct resolve_to_ublas<Transpose>
        { typedef boost::numeric::bindings::tag::row_major value; };
        
        template<class Tag>
        struct resolve_to_char;
        
        template<> struct resolve_to_char<NoTranspose>
        { static char value() { return 'N'; } };
        
        template<> struct resolve_to_char<Transpose>
        { static char value() { return 'T'; } };
        
#define MAKE_MY_GEMM(F, type) \
inline void mygemm(char const & c1, char const & c2, int const & m, int const & n, int const & k, \
type const * A, int const & lda, type const * B, int const & ldb, type * C, int const & ldc) \
{ \
    type one = type(1); \
    type zero = type(0); \
    F(&c1, &c2, &m, &n, &k, &one, A, &lda, B, &ldb, &zero, C, &ldc);\
}

    MAKE_MY_GEMM(dgemm_, double)
    MAKE_MY_GEMM(sgemm_, float)
    MAKE_MY_GEMM(cgemm_, std::complex<float>)
    MAKE_MY_GEMM(zgemm_, std::complex<double>)
    
#undef MAKE_MY_GEMM

    }
    
    template<typename T, typename MemoryBlock, class Tag1, class Tag2>
    void gemm(dense_matrix<T, MemoryBlock> const & m1, Tag1,
              dense_matrix<T, MemoryBlock> const & m2, Tag2,
              dense_matrix<T, MemoryBlock> & m3)
    {
        assert(detail::dims<Tag1>::second(m1) == detail::dims<Tag2>::first(m2));
        assert(detail::dims<Tag1>::first(m1) == m3.num_rows());
        assert(detail::dims<Tag2>::second(m2) == m3.num_cols());
        
        // This should work, but it doesn't. ;-)
//        boost::numeric::bindings::blas::gemm(typename dense_matrix<T,MemoryBlock>::value_type(1),
//                                             detail::TransConjWrapper<T, MemoryBlock, Tag1>(m1),
//                                             detail::TransConjWrapper<T, MemoryBlock, Tag2>(m2),
//                                             typename dense_matrix<T,MemoryBlock>::value_type(0),
//                                             m3);
        
        // The pedestrian way
        // This shouldn't be here, but it works.
        detail::mygemm(detail::resolve_to_char<Tag1>::value(), detail::resolve_to_char<Tag2>::value(),
                       m3.num_rows(), m3.num_cols(), detail::dims<Tag1>::second(m1),
                       &m1(0,0), m1.stride2(),
                       &m2(0,0), m2.stride2(),
                       &m3(0,0), m3.stride2());
    }
}

// this is needed only for the non-pedestrian approach above
namespace boost { namespace numeric { namespace bindings { namespace detail {
    
    using ::blas::detail::TransConjWrapper;
    
    // this is copied from dense_matrix_adaptor.hpp, with only the modification of the data_order
    
    template <typename T, typename MemoryBlock, typename Tag, typename Id, typename Enable>
    struct adaptor< TransConjWrapper<T,MemoryBlock,Tag>, Id, Enable>
    {
        typedef typename copy_const< Id, T >::type value_type;
        typedef std::ptrdiff_t size_type;
        typedef std::ptrdiff_t difference_type;
        
        typedef mpl::map<
        mpl::pair< tag::value_type,      value_type >,
        mpl::pair< tag::entity,          tag::matrix >,
        mpl::pair< tag::size_type<1>,    size_type >,
        mpl::pair< tag::size_type<2>,    size_type >,
        mpl::pair< tag::data_structure,  tag::linear_array >,
        mpl::pair< tag::data_order,      typename ::blas::detail::resolve_to_ublas<Tag>::value >,
        mpl::pair< tag::data_side,       tag::upper >,
        mpl::pair< tag::stride_type<1>,  tag::contiguous >,
        mpl::pair< tag::stride_type<2>,  difference_type >
        > property_map;
        
        static size_type size1( const Id& id ) {
            return id.ptr->num_rows();
        }
        
        static size_type size2( const Id& id ) {
            return id.ptr->num_cols();
        }
        
        static value_type* begin_value( Id& id ) {
            return &(*id.ptr->column(0).first);
        }
        
        static value_type* end_value( Id& id ) {
            return &(*(id.ptr->column(id.ptr->num_cols()-1).second-1));
        }
        
        static difference_type stride1( const Id& id ) {
            return id.ptr->stride1();
        }
        
        static difference_type stride2( const Id& id ) {
            return id.ptr->stride2();
        }
        
    };
    
    
}}}}

#endif
