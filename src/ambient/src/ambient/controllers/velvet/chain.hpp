namespace ambient { namespace controllers { namespace velvet {

    inline void* chain::operator new (size_t size){
        return ambient::bulk_pool.get<chain>();
    }

    inline void chain::operator delete (void* ptr){
    }

    inline chain::chain(cfunctor* f){
        this->push_back(f);
    }

    inline void chain::push_back(cfunctor* f){
        this->content.push_back(f);
        for(std::vector<revision*>::iterator it = f->derivatives.begin(); it != f->derivatives.end(); ++it) this->derivatives.push_back(*it);
        for(std::list<revision*>::iterator it = f->dependencies.begin(); it != f->dependencies.end(); ++it) this->dependencies.push_back(*it);
    }

    inline void chain::execute(){
        std::vector<cfunctor*>::iterator i;
        i = this->content.begin(); while(i != this->content.end()) (*i++)->computation();

        std::vector<revision*>::iterator it = this->derivatives.begin(); 
        while(it != this->derivatives.end()) (*it++)->reset_generator();

        i = this->content.begin(); while(i != this->content.end()) delete *i++;
    }

    inline bool chain::ready(){
        for(std::vector<revision*>::iterator it = this->dependencies.begin(); it != this->dependencies.end(); ++it) 
            if((*it)->generator != NULL) return false;
        return true;
    }

    inline bool chain::constrains(cfunctor* op){
        size_t check = op->dependencies.size();

        for(std::vector<revision*>::iterator it = this->derivatives.begin(); it != this->derivatives.end(); ++it)
        for(std::list<revision*>::iterator itd = op->dependencies.begin(); itd != op->dependencies.end(); ){
            if(*it == *itd) op->dependencies.erase(itd++);
            else itd++;
        }

        if(check != op->dependencies.size()) return true;
        return false;
    }

} } }
