namespace ambient { namespace controllers { namespace velvet {

    inline cfunctor::cfunctor(){
        ambient::controller.push(this);
    }

    inline cfunctor::~cfunctor(){
    //    this->grp->idle();
    } 

} } }
