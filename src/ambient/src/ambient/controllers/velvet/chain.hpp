namespace ambient { namespace controllers { namespace velvet {

    inline void* chain::operator new (size_t size){
        return ambient::bulk_pool.get<sizeof(chain)>();
    }

    inline void chain::operator delete (void* ptr){
    }

    inline chain::chain(cfunctor* op){
        this->content.reserve(2);
        this->push_back(op);
    }

    inline size_t chain::size(){
        return this->content.size();
    }

    inline void chain::push_back(cfunctor* op){
        this->content.push_back(op);
        op->tag(this);
    }

    inline void chain::execute(){
        std::vector<cfunctor*>::iterator i = this->content.begin(); 
        while(i != this->content.end())
            (*i++)->invoke();
    }

    inline bool chain::ready(){
        std::vector<cfunctor*>::iterator i = this->content.begin(); 
        while(i != this->content.end()) if(!(*i++)->ready(this)) return false;
        return true;
    }

} } }
