#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

#ifdef AMBIENT_TRACE_SYNC
#include <execinfo.h>
#endif

namespace ambient {

    inline void persist(history* o){ 
        controller.persist(o); 
    }

    template<typename T>
    inline void destroy(T* o){ 
        controller.destroy(o); 
    }

    inline bool verbose(){ 
        return rank.verbose;   
    }

    inline bool uniform(){
        return controller.uniform;
    }

    inline void log(const char* msg){
        if(ambient::rank()) printf("%s\n", msg);
    }

    inline void sync(){ 
#ifdef AMBIENT_TRACE_SYNC
        void *array[10];
        size_t size = backtrace(array, 10);
        backtrace_symbols_fd(array, size, 2);
#endif
        controller.flush();
        controller.clear();  
        bulk.drop();
    }

    inline void fuse(const history* src, history* dst){ 
        dst->fuse(src);
    }

}

#endif
