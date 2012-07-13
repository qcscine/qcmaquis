#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

namespace ambient {
    inline void playout()                { ambient::controller.flush();                  }
    inline void conditional_playout()    { ambient::controller.conditional_flush();      }
    inline void set_num_threads(size_t n){ ambient::controller.set_num_threads(n);       }
    inline size_t get_num_threads()      { return ambient::controller.get_num_threads(); }
    inline void mute()                   { ambient::controller.muted = true;             }
    inline void unmute()                 { ambient::controller.muted = false;            }
    inline bool verbose()                { return (rank() ? false : true);               }
}

#endif
