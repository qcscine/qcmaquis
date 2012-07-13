namespace ambient { namespace controllers { namespace velvet {

        inline cfunctor::cfunctor()
        : credit(0)
        {
            ambient::controller.push(this);
        }

        inline cfunctor::~cfunctor(){
            this->grp->idle();
        } 

        inline void cfunctor::complete(){
            return ambient::controller.atomic_complete(this);
        }

} } }
