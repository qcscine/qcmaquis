#ifndef __ALPS_DENSE_MATRIX_ADAPTOR_HPP__
#define __ALPS_DENSE_MATRIX_ADAPTOR_HPP__

#include <boost/numeric/bindings/detail/adaptor.hpp>
#include <boost/numeric/bindings/detail/if_row_major.hpp>

namespace maquis { namespace types {
    template <typename T, typename MemoryBlock> 
    class dense_matrix;
} }

//
// An adaptor for the matrix to the boost::numeric::bindings
//

namespace boost { namespace numeric { namespace bindings { namespace detail {
    
    template <typename T, typename MemoryBlock, typename Id, typename Enable>
    struct adaptor<::maquis::types::dense_matrix<T,MemoryBlock>, Id, Enable>
    {
        typedef typename copy_const< Id, T >::type              value_type;
        // TODO: fix the types of size and stride -> currently it's a workaround, since std::size_t causes problems with boost::numeric::bindings
        //typedef typename ::maquis::types::dense_matrix<T,Alloc>::size_type         size_type;
        //typedef typename ::maquis::types::dense_matrix<T,Alloc>::difference_type   difference_type;
        typedef std::ptrdiff_t  size_type;
        typedef std::ptrdiff_t  difference_type;

        typedef mpl::map<
            mpl::pair< tag::value_type,      value_type >,
            mpl::pair< tag::entity,          tag::matrix >,
            mpl::pair< tag::size_type<1>,    size_type >,
            mpl::pair< tag::size_type<2>,    size_type >,
            mpl::pair< tag::data_structure,  tag::linear_array >,
            mpl::pair< tag::data_order,      tag::column_major >,
            mpl::pair< tag::data_side,       tag::upper >,
            mpl::pair< tag::stride_type<1>,  tag::contiguous >,
            mpl::pair< tag::stride_type<2>,  difference_type >
        > property_map;

        static size_type size1( const Id& id ) {
            return id.num_rows();
        }

        static size_type size2( const Id& id ) {
            return id.num_cols();
        }

        static value_type* begin_value( Id& id ) {
            return &(*id.column(0).first);
        }

        static value_type* end_value( Id& id ) {
            return &(*(id.column(id.num_cols()-1).second-1));
        }

        static difference_type stride1( const Id& id ) {
            return id.stride1();
        }

        static difference_type stride2( const Id& id ) {
           return id.stride2();
        }

    };

}}}}

#endif //__ALPS_DENSE_MATRIX_ADAPTOR_HPP__
