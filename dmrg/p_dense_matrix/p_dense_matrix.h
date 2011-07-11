#ifndef __ALPS_P_DENSE_MATRIX_HPP__
#define __ALPS_P_DENSE_MATRIX_HPP__

#include "ambient/ambient.hpp"
#include "utils/zout.hpp"
#include "utils/matrix_vector_traits.h"

#include "p_dense_matrix/p_diagonal_matrix.h"

#include <ostream>

namespace blas {

    template <typename T, ambient::policy P>
    class p_dense_matrix : public ambient::livelong<p_dense_matrix<T,P>, P>
    {
    public:
        typedef p_dense_matrix<T, ambient::REPLICA> replica;
        typedef ambient::livelong<p_dense_matrix<T,P>, P> livelong;
        friend class p_dense_matrix<T, ambient::ANY>;
        friend class p_dense_matrix<T, ambient::REPLICA>;
        typedef T         value_type;                     // The type T of the elements of the matrix
        typedef size_t    size_type;                      // Unsigned integer type that represents the dimensions of the matrix
        typedef ptrdiff_t difference_type;                // Signed integer type to represent the distance of two elements in the memory

       ~p_dense_matrix();
        p_dense_matrix();
        p_dense_matrix(ambient::void_pt* p);              // proxy object construction
        p_dense_matrix(size_type rows, size_type cols, T init_v);
        p_dense_matrix(p_dense_matrix const& m);
  
        void swap(p_dense_matrix & r);
        friend void swap(p_dense_matrix & x, p_dense_matrix & y){ x.swap(y); }

        inline bool empty() const;
        inline value_type get_init_v() const;
        inline size_type num_rows() const;
        inline size_type num_cols() const;
        void resize(size_type rows, size_type cols);
        void remove_rows(size_type i, size_type k);
        void remove_cols(size_type j, size_type k);
        void nullcut();
        void clear();
        void touch() const;
       
        void inplace_conjugate();
        static p_dense_matrix<T,P> identity_matrix(size_type size);

        value_type& get(size_type i, size_type j) const;

        template<ambient::policy PR>                      operator p_dense_matrix<T,PR> ();
        template<ambient::policy PR> p_dense_matrix<T,P>& operator = (p_dense_matrix<T,PR>& rhs); // not ref. in orig
        p_dense_matrix<T,P>&                              operator = (const p_dense_matrix<T>& rhs);
        inline value_type&                                operator()(size_type i, size_type j);
        inline const value_type&                          operator()(size_type i, size_type j) const;
        p_dense_matrix<T,P>&                              operator += (const p_dense_matrix& rhs); 
        p_dense_matrix<T,P>&                              operator -= (const p_dense_matrix& rhs);
        template <typename T2> p_dense_matrix<T,P>&       operator *= (const T2& t);
        template <typename T2> p_dense_matrix<T,P>&       operator /= (const T2& t);

        value_type init_v; // value used for initialization
        
#ifdef HAVE_ALPS_HDF5
        void serialize(alps::hdf5::iarchive & ar) { std::cerr << "I don't do much." << std::endl; }
	void serialize(alps::hdf5::oarchive & ar) const { std::cerr << "I don't do much either." << std::endl; }
#endif        
    private:
        size_type rows;
        size_type cols;
    };

    template<typename T, ambient::policy P>
    struct associated_diagonal_matrix< p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };
    
    template<typename T, ambient::policy P>
    struct associated_real_diagonal_matrix< p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<typename detail::real_type<T>::type> type;
    };

    template<typename T, ambient::policy P>
    struct associated_vector<p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<T> type;
    };
    
    template<typename T, ambient::policy P>
    struct associated_real_vector<p_dense_matrix<T,P> >
    {
        typedef p_diagonal_matrix<typename detail::real_type<T>::type> type;
    };
} // namespace blas

#include "p_dense_matrix/p_dense_matrix.hpp"
#endif //__ALPS_DENSE_MATRIX_HPP__
