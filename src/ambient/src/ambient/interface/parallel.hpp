#ifndef AMBIENT_INTERFACE_PARALLEL
#define AMBIENT_INTERFACE_PARALLEL

namespace ambient{ 

    using ambient::models::velvet::history;
    using ambient::controllers::velvet::iteratable;
   
    template<typename T> class copy;
    template<typename T> class copy_atomic;

    template <typename T>
    class parallel : public iteratable<history>
    { 
    protected:
        typedef typename boost::intrusive_ptr<T> ptr;
        typedef typename ambient::info<T>::value_type value_type;
        friend void intrusive_ptr_add_ref<>(T* p);
        friend void intrusive_ptr_release<>(T* p);
        long references;

        inline parallel()
        : iteratable<history>(sizeof(value_type)), references(0) 
        {
        }

        inline parallel(const T& o)
        : iteratable<history>(sizeof(value_type)), references(1) // avoiding auto-deallocation
        {
            ambient::model.set_current_dim(this, ambient::model.get_current_dim(&o));
            if(!o.pt_atomic()) ambient::push< ambient::copy<value_type> >(*(T*)this, *(const T*)&o);
            else ambient::push< ambient::copy_atomic<value_type> >(*(T*)this, *(const T*)&o);
            this->references--; // push is done, can revert it back
        }

        inline void pt_set_dim(size_t x, size_t y){
            ambient::model.set_current_dim(this, dim2(x, y));
        }

        inline dim2 pt_mem_dim() const {
            return this->current->content->get_mem_dim();
        }

    public:
        inline bool pt_atomic() const {
            return ambient::model.is_atomic(this);
        }

        inline value_type& pt_fetch(size_t blocked_i, size_t blocked_j, 
                                    size_t element_i, size_t element_j){
            ambient::playout();
            return ((value_type*)ambient::controller.ufetch_block(*this->current, blocked_j, blocked_i))
                   [ element_j*this->pt_mem_dim().y+element_i ];
        }
    };

}

#endif
