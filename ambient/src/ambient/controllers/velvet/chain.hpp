namespace ambient { namespace controllers { namespace velvet {

    inline void* chain::operator new (size_t size){
        return ambient::bulk_pool.get<sizeof(chain)>();
    }

    inline void chain::operator delete (void* ptr){
    }

    inline chain::chain(cfunctor* op){
        this->push_back(op);
    }

    inline void chain::push_back(cfunctor* op){
        this->content.push_back(op);
        op->tag(this);
    }

    inline void chain::execute(){
        std::vector<cfunctor*>::iterator i = this->content.begin(); 
        while(i != this->content.end()){
            (*i)->computation();
            delete *i++;
        }
    }

    inline bool chain::ready(){
        std::vector<cfunctor*>::iterator i = this->content.begin(); 
        while(i != this->content.end()) if(!(*i++)->ready(this)) return false;
        return true;
    }

    inline bool chain::constrains(cfunctor* op){
        return op->match(this);
    }

} } }
