#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

namespace ambient {
    template<typename T>
    inline void destroy(T* o)                               { ambient::controller.destroy(o);          }
    inline bool verbose()                                   { return (rank() ? false : true);          }
    inline void sync(){ 
        ambient::controller.flush();
        ambient::controller.garbage.clear();
    } 
}

#endif
