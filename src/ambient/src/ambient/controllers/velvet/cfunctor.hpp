namespace ambient { namespace controllers { namespace velvet {

    class context_cilk{
    public:


    };    


    inline cfunctor::cfunctor()
    {
        this->credit = 0;
        ambient::controller.push(this);
        this->ctxt = new context_cilk();
    }

    inline cfunctor::~cfunctor(){
        this->grp->idle();
    } 

} } }
