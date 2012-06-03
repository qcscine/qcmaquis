namespace ambient { namespace models { namespace velvet {

    inline revision::revision(layout* l)
    : content(l), generator(NULL)
    {
    }

    inline revision::~revision(){
        delete this->content;
    }

    inline void revision::set_dim(dim2 dim){
        this->content->set_dim(dim);
    }

    inline dim2 revision::get_dim(){
        return this->content->get_dim();
    }

    inline size_t revision::id(){
        return this->content->id();
    }

    inline layout::entry* revision::block(size_t x, size_t y){
        return this->content->get(x, y);
    }

    inline void revision::add_modifier(sfunctor* m){
        this->modifiers.push_back(m);
    }

    inline std::list<sfunctor*>& revision::get_modifiers(){
        return this->modifiers;
    }

    inline group* revision::get_placement(){
        return this->content->placement;
    }

    inline void revision::set_placement(group* grp){
        this->content->placement = grp;
    }

    inline void revision::set_generator(sfunctor* m){
        this->generator = m;
    }

    inline sfunctor* revision::get_generator(){
        return this->generator;
    }

} } }
