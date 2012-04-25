#ifndef AMBIENT_INTERFACE_PARALLEL_T
#define AMBIENT_INTERFACE_PARALLEL_T

namespace ambient{ 

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
        : references(0)
        {
            this->t_size = sizeof(value_type);
        }

        parallel_t(const T& o)
        : references(0)
        {
            this->t_size = sizeof(value_type);
            this->pt_set_dim(o.get_dim().x, o.get_dim().y);
            ambient::push(ambient::copy_l<T>, ambient::copy_c<T>, 
                          *(T*)this, *(const T*)&o);
        }

        dim2 pt_mem_dim() const {
            return this->revision(0).get_layout().get_mem_dim();
        }

        void pt_set_dim(size_t x, size_t y = 1){
            this->set_dim(dim2(x, y));
        }

    public:
        template<typename O>
        void pt_set_init(void(*fp)(O&)){
            this->revision(0).init((void(*)())fp);
        }

        value_type& pt_fetch(size_t blocked_i, size_t blocked_j, 
                             size_t element_i, size_t element_j){
            ambient::playout();
            return ((value_type*)ambient::controller.ufetch_block(this->revision(0), blocked_i, blocked_j))
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
//    template <typename T>
//    models::imodel::revision& future(T& obj, int n = 2){
//        return ((models::imodel::object*)&obj)->revision(n);
//    }
    template <char R, typename T>
    models::imodel::reduction& reduced(T& obj){
        if(updated(obj).get_reduction() == NULL){
            if(R == '+') updated(obj).set_reduction();
        }
        return *updated(obj).get_reduction();
    }
    // }}}
}

#endif
