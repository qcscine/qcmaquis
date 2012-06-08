#ifndef AMBIENT_MODELS_VELVET_REDUCTION
#define AMBIENT_MODELS_VELVET_REDUCTION

namespace ambient { namespace models { namespace velvet {

    class reduction
    {
    public:
        class reductionq {
        public:
            reductionq();
            void push(layout::entry*);
        };
        reduction(revision*);
        ~reduction();
        layout::entry* block(size_t x, size_t y);
        layout::entry& operator()(size_t x, size_t y);
        std::vector< std::vector<reductionq*> > entries;
        revision* target;
    };

} } }

#endif
