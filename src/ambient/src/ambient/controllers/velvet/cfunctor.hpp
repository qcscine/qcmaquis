namespace ambient { namespace controllers { namespace velvet {

        inline cfunctor::cfunctor()
        : workload(1), credit(0)
        {
            pthread_mutex_init(&mutex, NULL);
            ambient::controller.push(this);
        }

        inline void cfunctor::check_complete(){
            lock(); 
            if(--workload == 0) return ambient::controller.atomic_complete(this);
            unlock();
        }

} } }
