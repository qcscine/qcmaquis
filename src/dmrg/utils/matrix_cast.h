#ifndef MATRIX_CAST_H
#define MATRIX_CAST_H

namespace blas{

    template<typename T, typename Memoryblock>
    class dense_matrix;

    template<typename T, enum Policy>
    class  p_dense_matrix;

    template<typename T>
    void cast(p_dense_matrix<T>& pm, dense_matrix<T> const& m)
    {
        const std::vector<typename dense_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows = m.num_rows();
        size_t num_cols = m.num_cols();
        ambient::push(ambient::cast_to_p_dense_l, ambient::cast_to_p_dense_c,v_ptr ,pm, num_rows, num_cols);
        ambient::playout();
    }

    template<typename T>
    void cast(dense_matrix<T>& m, p_dense_matrix<T> const& pm)
    {
        std::vector<typename dense_matrix<T>::value_type>* v_ptr = &m.get_values();
        size_t num_rows = m.num_rows();
        size_t num_cols = m.num_cols();
        ambient::push(ambient::cast_to_dense_l, ambient::cast_to_dense_c,v_ptr ,pm, num_rows, num_cols);
        ambient::playout();
    }

    template<typename T, typename V>
    T matrix_cast(V const& m_in)
    {
        T m_out(m_in.num_rows(), m_in.num_cols());    
        cast<typename V::value_type>(m_out,m_in);        
        return m_out;    
    }
}

#endif
