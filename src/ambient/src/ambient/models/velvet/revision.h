#ifndef AMBIENT_MODELS_VELVET_REVISION
#define AMBIENT_MODELS_VELVET_REVISION

namespace ambient { namespace models { namespace velvet {

    class revision
    {
    public:
        inline layout::entry& operator()(size_t x, size_t y){
            return this->block(x, y);
        } 
        inline layout& get_layout(){ return *content; }
        inline revision* get_parent(){ return parent; }
        inline revision(layout*);
        inline ~revision();
        inline layout::entry& block(size_t x, size_t y);
        inline void add_modifier(sfunctor* m);
        inline std::list<sfunctor*>& get_modifiers();
        inline size_t id();
        inline group* get_placement();
        inline void set_placement(group*);
        inline sfunctor* get_generator();
        inline void set_generator(sfunctor*);
        layout* const content;
        sfunctor* generator;
        std::list<sfunctor*> modifiers;
        revision* parent;
    };

} } }

#endif
