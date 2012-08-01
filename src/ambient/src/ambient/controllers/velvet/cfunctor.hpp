namespace ambient { namespace controllers { namespace velvet {

    inline void* cfunctor::operator new (size_t size){ 
        return ambient::bulk_pool.get<cfunctor>(); 
    }

    inline void cfunctor::operator delete (void* ptr){ }

    inline cfunctor::cfunctor(){
        ambient::controller.push(this);
    }

    inline cfunctor::~cfunctor(){
    //    this->grp->idle();
    } 

    inline bool cfunctor::ready(){
        std::list<revision*>::iterator it = this->dependencies.begin(); 
        while(it != this->dependencies.end()) if((*it++)->generator != NULL) return false;
        return true;
    }

} } }
