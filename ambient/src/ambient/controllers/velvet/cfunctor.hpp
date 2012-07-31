namespace ambient { namespace controllers { namespace velvet {

    inline void* cfunctor::operator new (size_t size){ 
        return ambient::bulk_pool.get<cfunctor>(); 
    }

    inline void cfunctor::operator delete (void* ptr){ }

    inline cfunctor::cfunctor()
    {
        this->credit = 0;
        ambient::controller.push(this);
    }

    inline cfunctor::~cfunctor(){
        this->grp->idle();
    } 

} } }
