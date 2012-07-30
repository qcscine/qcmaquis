namespace ambient { namespace controllers { namespace velvet {


    inline cfunctor::cfunctor()
    {
        this->credit = 0;
        ambient::controller.push(this);
    }

    inline cfunctor::~cfunctor(){
        this->grp->idle();
    } 

} } }
