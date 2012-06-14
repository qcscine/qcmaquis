#ifndef AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#define AMBIENT_CONTROLLERS_VELVET_CFUNCTOR
#include "ambient/utils/dim2.h"

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor : public models::velvet::sfunctor {
    public:
        virtual void weight()      = 0;
        virtual void logistics()   = 0;
        virtual void computation() = 0;
        virtual bool pretend()     = 0;
        virtual ~cfunctor(){
#ifdef AMBIENT_THREADS
            pthread_mutex_unlock(&mutex);
            pthread_mutex_destroy(&mutex);
#endif
        } 
        inline cfunctor();
        inline void check_complete();
#ifdef AMBIENT_THREADS
        inline void lock()                 { pthread_mutex_lock(&mutex);    }
        inline void unlock()               { pthread_mutex_unlock(&mutex);  }
#else
        inline void lock()                 { }
        inline void unlock()               { }
#endif
        inline void add_condition()        { lock(); workload++; unlock();  }
        inline void add_condition(size_t s){ lock(); workload+=s; unlock(); }
        inline size_t get_weight()         { return credit;                 }
        inline void set_weight(size_t c)   { credit = c;                    }
        long int workload; // signed for thread-safety
        size_t credit;
#ifdef AMBIENT_THREADS
        pthread_mutex_t mutex;
#endif
    };

} } }

#endif
