#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

// revision tracking mechanism (target selector)
namespace ambient { namespace models { namespace velvet {

    class history {
    public:
        inline void* operator new (size_t);
        inline void operator delete (void*);
        inline history(dim2,size_t);
        inline ~history();
        inline void clean();
        inline void add_state(revision* r);
        inline revision* back() const;
        inline size_t time() const;
        inline bool weak() const;
        std::vector<revision*> content;
        revision* current;
        size_t start;
        memspec spec;
        friend void intrusive_ptr_add_ref(history* p){ ++(p->references); } 
        friend void intrusive_ptr_release(history* p){ if(--(p->references) == 0) delete p; } 
    private:
        long references;
    };

} } }

#endif
