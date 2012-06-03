namespace ambient { namespace controllers { namespace velvet {

        inline void cfunctor::check_complete(){
            lock(); 
                if(--workload == 0) ::ambient::controller.atomic_complete(); 
            unlock();
        }

} } }
