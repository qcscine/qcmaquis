namespace ambient { namespace controllers { namespace velvet {

        inline cfunctor::cfunctor()
        {
            this->credit = 0;
            this->affinity = -1; 
            ambient::controller.push(this);
        }

        inline cfunctor::~cfunctor(){
            this->grp->idle();
        } 

} } }
