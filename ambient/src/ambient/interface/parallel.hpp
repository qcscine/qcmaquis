#ifndef AMBIENT_INTERFACE_PARALLEL
#define AMBIENT_INTERFACE_PARALLEL

namespace ambient{ 

    using ambient::models::velvet::history;
    using ambient::controllers::velvet::iteratable;
   
    template<typename T> class copy;

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
            if(o.get_dim().y == 0 || o.get_dim().x == 0) printf("Evil is here! is copy!\n");
            this->pt_set_dim(o.get_dim().x, o.get_dim().y);
            ambient::push< ambient::copy<value_type> >(*(T*)this, *(const T*)&o);
            this->references--; // push is done, can revert it back
        }

        inline dim2 pt_mem_dim() const {
            return this->current->content->get_mem_dim();
        }

        inline dim2 pt_grid_dim() const {
            return this->current->content->grid_dim;
        }

        inline void pt_set_dim(size_t x, size_t y){
            if(x == 0 || y == 0) printf("Evil is outside!\n");
            this->set_dim(dim2(x, y));
        }

    public:
        inline value_type& pt_fetch(size_t blocked_i, size_t blocked_j, 
                                    size_t element_i, size_t element_j){
            ambient::playout();
            return ((value_type*)ambient::controller.ufetch_block(*this->current, blocked_j, blocked_i))
                   [ element_j*this->pt_mem_dim().y+element_i ];
        }
    };

}

#endif
