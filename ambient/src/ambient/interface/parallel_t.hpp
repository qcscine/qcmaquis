// nested inside ambient.hpp in ambient namespace
#define BOOST_SP_NO_SP_CONVERTIBLE
#include <boost/intrusive_ptr.hpp>

#define pinned ambient::models::ambient_pin* ,

template <typename T>
void copy_l(maquis::types::p_dense_matrix_impl<T>& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a);

template <typename T>
void copy_c(maquis::types::p_dense_matrix_impl<T>& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a);

template <typename T>
class parallel_t : public models::v_model::object 
{ 
protected:
    typedef typename boost::intrusive_ptr<T> ptr;
    typedef typename ambient::models::info<T>::value_type value_type;
    friend void intrusive_ptr_add_ref<>(T* p);
    friend void intrusive_ptr_release<>(T* p);
    long references;

    parallel_t()
    :references(0)
    {
        this->t_size = sizeof(value_type);
    }

    parallel_t(const T& o)
    :references(0)
    {
        this->t_size = sizeof(value_type);
        ambient::push(ambient::copy_l<value_type>, ambient::copy_c<value_type>, 
                      *(T*)this, *(const T*)&o);
    }

    dim2 pt_mem_dim() const {
        return this->revision(0).get_layout().get_mem_dim();
    }

    void pt_set_dim(size_t x, size_t y = 1){
        this->dim.x = x;
        this->dim.y = y;
    }

public:
    template<typename O>
    void pt_set_init(void(*fp)(O&)){
        //this->init(new models::operation(fp, (T*)this));
    }

    value_type& pt_fetch(size_t blocked_i, size_t blocked_j, 
                         size_t element_i, size_t element_j){
        ambient::playout();
        return ((value_type*)ambient::controller.fetch_block(this->revision(0), blocked_i, blocked_j))
               [ element_j*this->pt_mem_dim().y+element_i ];
    }
};

// {{{ aliases for use inside kernels
template <typename T>
models::imodel::revision& current(T& obj){
    return ((models::imodel::object*)&obj)->revision(0);
}
template <typename T>
models::imodel::revision& updated(T& obj){
    return ((models::imodel::object*)&obj)->revision(1);
}
template <typename T>
models::imodel::revision& future(T& obj, int n = 2){
    return ((models::imodel::object*)&obj)->revision(n);
}
template <char R, typename T>
models::imodel::revision& reduced(T& obj){
    //if(future(obj).reduction == NULL){ // 1 reduction at a time
        //if(R == '+') future(obj).reduce(plus_reduce<T>);
    //}
    return future(obj);
}
// }}}
