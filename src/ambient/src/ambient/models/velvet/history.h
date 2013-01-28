#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

// revision tracking mechanism (target selector)
namespace ambient { namespace models { namespace velvet {

    class history {
    public:
        void* operator new (size_t);
        void operator delete (void*);
        history(dim2,size_t);
        void add_state(void* g);
        void fuse(const history* src);
        revision* back() const;
        size_t time() const;
        bool weak() const;
        std::vector<revision*> content;
        revision* current;
        memspec spec;
        size_t clock;
    };

} } }

#endif
