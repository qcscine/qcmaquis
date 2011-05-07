#include "p_dense_matrix/p_dense_matrix.h"
#define STRONG_BARRIER MPI_Barrier(MPI_COMM_WORLD);fflush(stdout);


namespace blas {

    #define self this->self
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>::~p_dense_matrix(){ // #destructor
    }
    template <typename T, ambient::policy P> // proxy object construction
    p_dense_matrix<T,P>::p_dense_matrix(ambient::void_pt* p): livelong(p){ 
    }
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>::p_dense_matrix(size_type rows = 0, size_type cols = 0, T init_value = T() )
    : livelong(rows, cols, init_value){
        self->rows = rows;
        self->cols = cols;
        this->set_breakdown();
    }

    template <typename T, ambient::policy P>
    template <typename TM, ambient::policy PM>
    p_dense_matrix<T,P>::p_dense_matrix(p_dense_matrix<TM,PM> const& m)
    : livelong(m){
        self->rows = m.num_rows();
        self->cols = m.num_cols();
        this->set_breakdown();
        ambient::push(ambient::copy_l_kernel, ambient::copy_c_kernel, *this, m);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::swap(p_dense_matrix & r)
    {
        assert(false);
        std::swap(self->breakdown(), r.thyself->breakdown());
        std::swap(self->rows,        r.thyself->rows);
        std::swap(self->cols,        r.thyself->cols);
    }

    template <typename T, ambient::policy P>
    inline bool p_dense_matrix<T,P>::empty() const { return (self->rows == 0 || self->cols == 0); }

    template <typename T, ambient::policy P>
    inline size_t& p_dense_matrix<T,P>::num_rows() const { return self->rows; }

    template <typename T, ambient::policy P>
    inline size_t& p_dense_matrix<T,P>::num_cols() const { return self->cols; }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::nullcut(){
        if(!self->breakdown()->is_loose())
            ambient::push(ambient::nullcut_l_kernel, ambient::nullcut_c_kernel, *this, self->rows, self->cols);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::clear(){
        self->rows = self->cols = 0;
        this->nullcut();
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::resize(size_type rows, size_type cols)
    {
        self->rows = rows; 
        self->cols = cols;
        if(rows <= 0){ printf("%d %d\n", rows, cols); void* a = malloc(1); free(a); free(a); assert(rows > 0); }
        assert(cols > 0);
        if(self->breakdown()->is_loose()) self->breakdown()->set_dim(ambient::dim2(cols,rows));
        else ambient::push(ambient::resize_l_kernel, ambient::resize_c_kernel, *this, self->rows, self->cols);
        this->nullcut();
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::remove_rows(size_type i, size_type k = 1)
    {
        assert( i+k <= self->rows );
        if(self->breakdown()->is_loose()) return this->resize(self->rows-k, self->cols);
        ambient::push(ambient::remove_rows_l_kernel, ambient::remove_rows_c_kernel, *this, i, k);
        this->resize(self->rows - k, self->cols);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::remove_cols(size_type j, size_type k = 1)
    {
        assert( j+k <= self->cols );
        if(self->breakdown()->is_loose()) return this->resize(self->rows, self->cols-k);
        ambient::push(ambient::remove_cols_l_kernel, ambient::remove_cols_c_kernel, *this, j, k);
        this->resize(self->rows, self->cols - k);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::inplace_conjugate()
    {
        assert(false);
    }

    /*template <typename T> // alternative version // to use in future
    inline T& p_dense_matrix<T>::get(size_type i, size_type j) const
    {
        assert(i < self->rows);
        assert(j < self->cols);
        ambient::playout();
        int block_i = i / (self->breakdown()->get_mem_t_dim().y);
        int block_j = j / (self->breakdown()->get_mem_t_dim().x);
        int element_i = i % (self->breakdown()->get_mem_t_dim().y);
        int element_j = j % (self->breakdown()->get_mem_t_dim().x);
        if(self->breakdown()->block(block_i,block_j)->available())
            return *(T*)(*self->breakdown())(block_i, block_j).element(element_i, element_j);
        else //return *(new T()); //using default value of T
            throw ambient::core::remote_memory_e();
    }*/

    template <typename T, ambient::policy P>
    inline T& p_dense_matrix<T,P>::get(size_type i, size_type j) const
    {
        static double zero = 13;
        if(i >= self->rows){ printf("Accesing %d %d (of %d %d)\n", (int)i, (int)j, (int)self->rows, (int)self->cols); void* a = malloc(1); free(a); free(a); assert(i < self->rows); }
        assert(j < self->cols);
        if(ambient::occupied()) return zero; // for stream reader
        if(self->breakdown()->loose){
//            printf("Warning: accessing loose data...\n");
            ambient::push(ambient::nullcut_l_kernel, ambient::nullcut_c_kernel, *this, self->rows, self->cols);
        }
        ambient::playout();
        int block_i = i / (self->breakdown()->get_mem_t_dim().y);
        int block_j = j / (self->breakdown()->get_mem_t_dim().x);
        int element_i = i % (self->breakdown()->get_mem_t_dim().y);
        int element_j = j % (self->breakdown()->get_mem_t_dim().x);
        T& value = *(T*)(*self->breakdown())(block_i, block_j).element(element_i, element_j);
        ambient::world_loop();
        return value;
    }

    template <typename T, ambient::policy P>
    inline T&            p_dense_matrix<T,P>::operator () (size_type i, size_type j){ return this->get(i,j); }
    template <typename T, ambient::policy P>
    inline const T&      p_dense_matrix<T,P>::operator () (size_type i, size_type j) const { return this->get(i,j); }
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator += (const p_dense_matrix& rhs){ return(*this = *this + rhs); }
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator -= (const p_dense_matrix& rhs){ return (*this = *this - rhs); }
    template <typename T, ambient::policy P>
    template <typename T2>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator *= (const T2& t){ return (*this = *this * t); }
    template <typename T, ambient::policy P>
    template <typename T2>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator /= (const T2& t){ return (*this = *this * (1/t)); }
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator  = (const p_dense_matrix<T>& rhs){
        if(rhs.thyself->breakdown()->is_proxy()){
            ambient::pin(rhs, *this); // no copying - pinning profile
        }else{
            if(this->num_rows() != rhs.num_rows() || this->num_cols() != rhs.num_cols()) this->resize(rhs.num_rows(), rhs.num_cols());
            ambient::push(ambient::copy_l_kernel, ambient::copy_c_kernel, *this, rhs);
        }
        return *this;
    }

    template <typename T, ambient::policy P>
    template <ambient::policy PR>
    p_dense_matrix<T,P>::operator p_dense_matrix<T,PR>(){
        return (*(p_dense_matrix<T,PR>*)this);
    }
    template <typename T, ambient::policy P>
    template <ambient::policy PR>
    p_dense_matrix<T,P>& p_dense_matrix<T,P>::operator = (p_dense_matrix<T,PR>& rhs){
        return (*this = *(const p_dense_matrix<T>*)&rhs);
    }
    template <typename T, ambient::policy P>
    const p_dense_matrix<T>& operator + (const p_dense_matrix<T,P>& a, const p_dense_matrix<T,P>& b){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::mem_bound_l_kernel, ambient::add_c_kernel, a, b); 
        breakdown(out).set_init(ambient::nullify<T>);
        return out; 
    }
    template <typename T, ambient::policy P>
    const p_dense_matrix<T>& operator - (const p_dense_matrix<T,P>& a, const p_dense_matrix<T,P>& b){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::mem_bound_l_kernel, ambient::sub_c_kernel, a, b); 
        breakdown(out).set_init(ambient::nullify<T>);
        return out; 
    }
    template<typename T, ambient::policy P>
    const p_dense_matrix<T>& operator * (const p_dense_matrix<T,P>& lhs, const p_dense_matrix<T,P>& rhs){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::gemm_l_kernel, ambient::gemm_c_kernel, lhs, rhs);
        breakdown(out).set_init(ambient::nullify<T>);
        return out; 
    }
    template<typename T, ambient::policy P, typename T2>
    const p_dense_matrix<T>& operator * (const p_dense_matrix<T,P>& m, const T2& t){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::scale_l_kernel, ambient::scale_c_kernel, m, t); 
        breakdown(out).set_init(ambient::nullify<T>);
        return out; 
    }
    template<typename T, ambient::policy P, typename T2>
    const p_dense_matrix<T>& operator * (const T2& t, const p_dense_matrix<T,P>& m){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::scale_l_kernel, ambient::scale_c_kernel, m, t); 
        breakdown(out).set_init(ambient::nullify<T>);
        return out; 
    }
//////////////////////////////////// AMBIENT PART ////////////////////////////////////////////////////

    /*template <typename T>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T> const& m)
    {
        ambient::playout();
        STRONG_BARRIER
        for(typename p_dense_matrix<T>::size_type i=0; i< m.num_rows(); ++i)
        {
            STRONG_BARRIER
            for(typename p_dense_matrix<T>::size_type j=0; j < m.num_cols(); ++j){
                STRONG_BARRIER
                try{
                    m(i,j); printf("%.2f ", m(i,j)); // faster than cout
                }catch(...){ usleep(10*1000); }      // 10 ms sleep
                STRONG_BARRIER
            }
            if(ambient::is_master()) printf("\n");
            STRONG_BARRIER
        }
        return o;
    }*/

    template <typename T, ambient::policy P>
    std::ostream& operator << (std::ostream& o, p_dense_matrix<T,P> const& m)
    {
        STRONG_BARRIER
        for(typename p_dense_matrix<T,P>::size_type i=0; i< m.num_rows(); ++i)
        {
            STRONG_BARRIER
            for(typename p_dense_matrix<T,P>::size_type j=0; j < m.num_cols(); ++j){
                STRONG_BARRIER
                if(ambient::is_master()) printf("%.2f	", m(i,j));
                else m(i,j); // just touch
                STRONG_BARRIER
            }
            if(ambient::is_master()) printf("\n");
            STRONG_BARRIER
        }
        return o;
    }
    #undef self
}
