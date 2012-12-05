#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

namespace ambient {
    template<typename T>
    inline void destroy(T* o){ controller.destroy(o); }
    inline bool verbose()    { return rank.verbose;   }
    inline void sync()       { controller.flush();
                               controller.clear();    } 
}

#endif
