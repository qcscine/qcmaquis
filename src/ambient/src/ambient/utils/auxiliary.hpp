#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

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

    inline void log(const char* msg){
        if(ambient::rank()) printf("%s\n", msg);
    }

    inline void sync(){ 
        controller.flush();
        controller.clear();  
        bulk.drop();
    }

    inline void fuse(const history* src, history* dst){ 
        dst->fuse(src);
    }

}

#endif
