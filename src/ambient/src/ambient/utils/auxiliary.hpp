#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

namespace ambient {
    template<typename T>
    inline void destroy(T* o)                               { ambient::controller.destroy(o);          }
    inline void conditional_sync()                          { ambient::controller.conditional_flush(); }
    inline void mute()                                      { ambient::controller.muted = true;        }
    inline void unmute()                                    { ambient::controller.muted = false;       }
    inline bool verbose()                                   { return (rank() ? false : true);          }
    inline void sync(){ 
        ambient::controller.schedule();             
        if(ambient::controller.chains.size() == 1){
            ambient::controller.execute(ambient::controller.chains.pick());
            ambient::controller.chains.reset();
        }else{
            ambient::controller.flush();
            ambient::controller.garbage.clear();
        }
    } 
}

#endif
