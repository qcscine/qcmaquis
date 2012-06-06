#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

namespace ambient { namespace models { namespace velvet {

    // revision tracking mechanism (target selector)
    class history {
    protected:
        inline history(size_t ts);
    public:
        inline ~history();
        inline revision& add_state(layout* l);
        inline size_t get_t_size() const; // is it used?
        inline size_t time() const;
        inline revision* back() const;
        inline dim2 get_cached_dim() const; // model
        inline void cache_dim(dim2); // model
        std::vector<revision*> content;
        size_t t_size;
        dim2   dim;
        revision* current;
    };

} } }

#endif
