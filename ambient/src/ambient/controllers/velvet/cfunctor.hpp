namespace ambient { namespace controllers { namespace velvet {

        inline cfunctor::cfunctor()
        : workload(1), credit(0)
        {
            this->derivative = NULL;
            pthread_mutex_init(&mutex, NULL);
            ambient::controller.push(this);
        }

        inline cfunctor::~cfunctor(){
            this->grp->idle();
            //if(this->derivative){
            //    this->derivative->reset_generator();
            //    ambient::controller.atomic_receive(*this->derivative, 0, 0);
            //}
            pthread_mutex_unlock(&mutex);
            pthread_mutex_destroy(&mutex);
        } 

        inline void cfunctor::check_complete(){
            lock(); 
            if(--workload == 0) return ambient::controller.atomic_complete(this);
            unlock();
        }

} } }
