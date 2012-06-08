#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

// revision tracking mechanism (target selector)
namespace ambient { namespace models { namespace velvet {

    class history {
    protected:
        inline history();
    public:
        inline ~history();
        inline void clean();
        inline revision& add_state(layout* l);
        inline size_t time() const;
        inline revision* back() const;
        inline dim2 get_cached_dim() const; // model
        inline void cache_dim(dim2); // model
        std::vector<revision*> content;
        revision* current;
        size_t start;
        dim2 dim;
    };

} } }

#endif
