#ifndef MATRIX_CAST_H
#define MATRIX_CAST_H

#include "types/dense_matrix/dense_matrix.h"
#include "types/dense_matrix/matrix_interface.hpp"
#include "types/dense_matrix/diagonal_matrix.h"
// C - This file must include into the shared memory version so ifdef
#ifdef MPI_PARALLEL 

#include "types/p_dense_matrix/p_dense_matrix.h"
#include "types/p_dense_matrix/p_diagonal_matrix.h"

namespace maquis
{
    namespace traits
    {
    template<typename T>
    void cast(maquis::types::p_dense_matrix<T>& pm, maquis::types::dense_matrix<T> const& m)
    {
        const std::vector<typename maquis::types::dense_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows = m.num_rows();
        size_t num_cols = m.num_cols();
        size_t lda = m.stride2();
        ambient::push(ambient::cast_to_p_dense_l<T>, ambient::cast_to_p_dense_c<T>, v_ptr, *pm.impl, num_rows, num_cols, lda);
        ambient::playout();
    }

    template<typename T>
    void cast(maquis::types::dense_matrix<T>& m, maquis::types::p_dense_matrix<T> const& pm)
    {
        std::vector<typename maquis::types::dense_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows = m.num_rows();
        size_t num_cols = m.num_cols();
        ambient::push(ambient::cast_to_dense_l<T>, ambient::cast_to_dense_c<T>, v_ptr, *pm.impl, num_rows, num_cols);
        ambient::playout();
    }

    template<typename T>
    void cast(maquis::types::p_diagonal_matrix<T>& pm, maquis::types::diagonal_matrix<T> const& m)
    {
        const std::vector<typename maquis::types::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows(m.num_rows());
        size_t num_cols(1);
        ambient::push(ambient::cast_to_p_dense_l<T>, ambient::cast_to_p_dense_c<T>, v_ptr, *pm.get_data().impl, num_rows, num_cols, num_rows);
        ambient::playout();
    }

    template<typename T>
    void cast(maquis::types::diagonal_matrix<T>& m, maquis::types::p_diagonal_matrix<T> const& pm)
    {
        std::vector<typename maquis::types::diagonal_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows = m.num_rows();
        size_t num_cols(1);
        ambient::push(ambient::cast_to_dense_l<T>, ambient::cast_to_dense_c<T>, v_ptr, *pm.get_data().impl, num_rows, num_cols);
        ambient::playout();
    }

    template<typename T>
    void cast(maquis::types::p_dense_matrix<T>& pm0, maquis::types::p_dense_matrix<T> const& pm1)
    {
        printf("You are casting to same types : p_Dense_Matrix \n");
    }

    template<typename T>
    void cast(maquis::types::p_diagonal_matrix<T>& pm0, maquis::types::p_diagonal_matrix<T> const& pm1)
    {
        printf("You are casting to same types : p_Diagonal_Matrix \n");
    }

    } // name space traits
} //end name space maquis
#endif

namespace maquis{
    namespace traits{

    template<typename T>
    void cast(maquis::types::dense_matrix<T>& m0, maquis::types::dense_matrix<T> const& m1)
    {
        printf("You are casting to same types : Dense_Matrix \n");
        m0 = m1; // bad
    }

    template<typename T>
    void cast(maquis::types::diagonal_matrix<T>& m0, maquis::types::diagonal_matrix<T> const& m1)
    {
        printf("You are casting to same types : Diagonal_Matrix \n");
        m0 = m1; // bad
    }

    template<typename T, typename V> // p_dense2dense or dense2p_dense
    T matrix_cast(V const& m_in)
    {
        T m_out(m_in.num_rows(), m_in.num_cols());    
        cast(m_out,m_in);        
        return m_out;    
    }

    } // name space traits
} //end name space maquis

#ifdef MPI_PARALLEL
    template<typename T>
    bool operator == (maquis::types::dense_matrix<T> const & m, maquis::types::p_dense_matrix<T> const & pm) 
    {
        int *n = new int(1); // C - default matrices are equals 
        maquis::types::p_dense_matrix<T> pma(maquis::traits::matrix_cast<maquis::types::p_dense_matrix<T> >(m));
        ambient::push(ambient::validation_l<T>, ambient::validation_c<T>, *pm.impl, *pma.impl, n); 
        ambient::playout();
        bool b(*n);
        delete n;            
        return b;
    }

    template<typename T>
    bool operator == (maquis::types::diagonal_matrix<T> const & m, maquis::types::p_diagonal_matrix<T> const & pm) 
    {
        int *n = new int(1); // C - default matrices are equals 
        maquis::types::p_diagonal_matrix<T>  pma(maquis::traits::matrix_cast<maquis::types::p_diagonal_matrix<T> >(m));
        ambient::push(ambient::validation_l<T>, ambient::validation_c<T>, *pm.get_data().impl, *pma.get_data().impl, n); 
        ambient::playout();
        bool b(*n);
        delete n;            
        return b;
    }

    template<typename T>
    bool operator == (maquis::types::p_diagonal_matrix<T> const & pm, maquis::types::diagonal_matrix<T> const & m) 
    {
        return (m == pm);
    }

    template<typename T>
    bool operator == (maquis::types::p_dense_matrix<T> const & pm, maquis::types::dense_matrix<T> const & m) 
    {
        return (m == pm);
    }
#endif

#endif
