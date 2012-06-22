#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

// revision tracking mechanism (target selector)
namespace ambient { namespace models { namespace velvet {

    class history {
    protected:
        inline history(dim2);
    public:
        inline ~history();
        inline void clean();
        inline revision& add_state(revision* r);
        inline size_t time() const;
        inline revision* back() const;
        inline bool weak() const;
        std::vector<revision*> content;
        revision* current;
        size_t start;
        memspec spec;
    };

} } }

#endif
