#ifndef AMBIENT_CONTROLLERS_VELVET_ITERATABLE
#define AMBIENT_CONTROLLERS_VELVET_ITERATABLE

extern pthread_key_t pthread_tid;
#define GET_TID 0 //*(size_t*)pthread_getspecific(pthread_tid)

namespace ambient { namespace controllers { namespace velvet {

    using ambient::models::velvet::revision;
    using ambient::models::velvet::layout;

    class fast_revision
    {
    public:
        inline layout::entry& operator()(size_t x, size_t y);
        inline layout& get_layout();
    };

    class slow_revision
    {
    public:
        inline layout::entry& operator()(size_t x, size_t y); // should be operator
        inline layout& get_layout();
    };

    template<class T>
    class iteratable : public T
    {
    protected:
        inline iteratable();
        inline ~iteratable();
    public:
        inline revision& state(size_t offset) const;
        inline revision& ui_l_revision_0() const {
            return *this->content[this->thread_revision_base[GET_TID]]; 
        }
        inline revision& ui_l_revision_1() const { // only for debug for now
            return *this->content[this->thread_revision_base[GET_TID] + 1]; 
        }
        inline fast_revision& ui_c_revision_0() const {
            return *(fast_revision*)this->content[this->thread_revision_base[GET_TID]]; 
        }
        inline slow_revision& ui_c_revision_1() const { 
            return *(slow_revision*)this->content[this->thread_revision_base[GET_TID] + 1]; 
        }
        inline size_t get_thread_revision_base() const;
        inline void set_thread_revision_base(size_t);

        size_t* thread_revision_base;
    };

} } }

#endif
