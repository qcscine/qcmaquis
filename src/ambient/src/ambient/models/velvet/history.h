#ifndef AMBIENT_MODELS_VELVET_HISTORY
#define AMBIENT_MODELS_VELVET_HISTORY

namespace ambient { namespace models { namespace velvet {

    // revision tracking mechanism (target selector)
    class history {
    protected:
        inline history(size_t ts);
    public:
        inline ~history();
        inline void add_state(layout* l);
        inline size_t get_t_size() const; // is it used?
        inline size_t time() const;
        inline dim2 get_dim() const;
        inline void set_dim(dim2);
        inline revision* back();
        std::vector<revision*> content;
        size_t t_size;
        revision* current;
    };

} } }

#endif
