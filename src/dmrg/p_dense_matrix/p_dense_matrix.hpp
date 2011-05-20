#include "p_dense_matrix/p_dense_matrix.h"
#define STRONG_BARRIER MPI_Barrier(MPI_COMM_WORLD);fflush(stdout);


namespace blas {

    #define self this->self
    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>::~p_dense_matrix(){ // #destructor
//        if(this->is_loose_copied() && !this->is_loose_copy() && !this->is_loose())
//            this->resize(self->duplicant->num_rows(), self->duplicant->num_cols());
    }
    template <typename T, ambient::policy P> // proxy object construction
    p_dense_matrix<T,P>::p_dense_matrix(ambient::void_pt* p): livelong(p){}

    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>::p_dense_matrix(size_type rows = 0, size_type cols = 0, T init_v = T() )
    : livelong(rows, cols, init_v){
        self->rows   = rows;
        self->cols   = cols;
        self->init_v = init_v;
        this->set_breakdown();
    }

    template <typename T, ambient::policy P>
    p_dense_matrix<T,P>::p_dense_matrix(const p_dense_matrix& m) // copy constructor (to be optimized)
    : livelong(m){
        self->rows   = m.num_rows();
        self->cols   = m.num_cols();
        self->init_v = m.init_v;
        this->set_breakdown();
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::swap(p_dense_matrix & r)
    {
        assert(false);
        std::swap(this->breakdown(), r.breakdown());
        std::swap(self->rows,        r.thyself->rows);
        std::swap(self->cols,        r.thyself->cols);
        std::swap(self->init_v,      r.thyself->init_v);
    }

    template <typename T, ambient::policy P>
    inline bool p_dense_matrix<T,P>::empty() const { return (self->rows == 0 || self->cols == 0); }

    template <typename T, ambient::policy P>
    inline size_t p_dense_matrix<T,P>::num_rows() const { return self->rows; }

    template <typename T, ambient::policy P>
    inline size_t p_dense_matrix<T,P>::num_cols() const { return self->cols; }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::nullcut(){
        ambient::push(ambient::nullcut_l, ambient::nullcut_c, *this, self->rows, self->cols);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::clear(){
        self->rows = self->cols = 0;
        if(!self->is_loose()) this->nullcut();
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::touch() const {
        self->abstract = false;
        ambient::push(ambient::touch_l, ambient::touch_c, *this);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::resize(size_type rows, size_type cols)
    {
        self->rows = rows; 
        self->cols = cols;
        assert(cols > 0);
        assert(rows > 0);
        if(!self->is_loose()){ 
            ambient::push(ambient::resize_l, ambient::resize_c, 
                          *this, self->rows, self->cols);
        }else
            self->breakdown()->set_dim(ambient::dim2((unsigned int)rows, (unsigned int)cols));
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::remove_rows(size_type i, size_type k = 1)
    {
        assert( i+k <= self->rows );
        if(self->is_loose()) return this->resize(self->rows-k, self->cols);
        ambient::push(ambient::remove_rows_l, ambient::remove_rows_c, *this, i, k);
        this->resize(self->rows - k, self->cols);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::remove_cols(size_type j, size_type k = 1)
    {
        assert( j+k <= self->cols );
        if(self->is_loose()) return this->resize(self->rows, self->cols-k);
        ambient::push(ambient::remove_cols_l, ambient::remove_cols_c, *this, j, k);
        this->resize(self->rows, self->cols - k);
    }

    template <typename T, ambient::policy P>
    void p_dense_matrix<T,P>::inplace_conjugate()
    {
        assert(false);
    }

    template<typename T, ambient::policy P>
    p_dense_matrix<T,P> p_dense_matrix<T,P>::identity_matrix(size_type size)
    {
        p_dense_matrix<T,P> ret(size, size);
        ret.set_init(ambient::identity_i<T>);
        return ret;
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
        assert(i < self->rows);
        assert(j < self->cols);

        if(ambient::access.write_only_marked()){
        // can write into some buffer
            if(self->modifier == NULL){
                self->modifiers.push(std::vector< std::pair<std::pair<size_t,size_t>,void*> >());
                ambient::push(ambient::apply_writes_l<T>, ambient::apply_writes_c<T>, *this);
                self->modifier = &self->modifiers.back();
            }
            void* value = calloc(1, sizeof(T));
            self->modifier->push_back(std::pair<std::pair<size_t,size_t>,void*>(std::pair<size_t,size_t>(i,j), value));
            return *(T*)value;
        }

        if(ambient::occupied()) assert(false);
        if(self->is_abstract()) this->touch();
        ambient::playout();
        int block_i = i / (this->breakdown()->get_mem_t_dim().y);
        int block_j = j / (this->breakdown()->get_mem_t_dim().x);
        int element_i = i % (this->breakdown()->get_mem_t_dim().y);
        int element_j = j % (this->breakdown()->get_mem_t_dim().x);
        T& value = *(T*)(*this->breakdown())(block_i, block_j).element(element_i, element_j);
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
        if(rhs.breakdown()->is_proxy()){
            ambient::pin(rhs, *this); // no copying - pinning profile
        }else{
            if(this->num_rows() != rhs.num_rows() || this->num_cols() != rhs.num_cols()) 
                this->resize(rhs.num_rows(), rhs.num_cols());
            ambient::push(ambient::copy_l, ambient::copy_c, *this, rhs);
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
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::mem_bound_l, ambient::add_c, a, b); 
        out.set_init(ambient::null_i<T>);
        return out; 
    } // returns non-ref in orig
    template <typename T, ambient::policy P>
    const p_dense_matrix<T>& operator - (const p_dense_matrix<T,P>& a, const p_dense_matrix<T,P>& b){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::mem_bound_l, ambient::sub_c, a, b); 
        out.set_init(ambient::null_i<T>);
        return out; 
    }
    template<typename T, ambient::policy P>
    const p_dense_matrix<T>& operator * (const p_dense_matrix<T,P>& lhs, const p_dense_matrix<T,P>& rhs){
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::gemm_l, ambient::gemm_c, lhs, rhs);
        out.set_init(ambient::null_i<T>);
        return out; 
    }
    template<typename T, ambient::policy P, typename T2>
    const p_dense_matrix<T>& operator * (const p_dense_matrix<T,P>& m, const T2& t){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::scale_l, ambient::scale_c, m, t); 
        out.set_init(ambient::null_i<T>);
        return out; 
    }
    template<typename T, ambient::policy P, typename T2>
    const p_dense_matrix<T>& operator * (const T2& t, const p_dense_matrix<T,P>& m){ 
        p_dense_matrix<T>& out = ambient::push< p_dense_matrix<T> >(ambient::scale_l, ambient::scale_c, m, t); 
        out.set_init(ambient::null_i<T>);
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
