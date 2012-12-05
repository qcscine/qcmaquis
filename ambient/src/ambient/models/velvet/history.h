#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

// revision tracking mechanism (target selector)
namespace ambient { namespace models { namespace velvet {

    class history {
    public:
        void* operator new (size_t);
        void operator delete (void*);
        history(dim2,size_t);
        ~history();
        void add_state(revision* r);
        revision* back() const;
        size_t time() const;
        bool weak() const;
        std::vector<revision*> content;
        revision* current;
        memspec spec;
    };

} } }

#endif
