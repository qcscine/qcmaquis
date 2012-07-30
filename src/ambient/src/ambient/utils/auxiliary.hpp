#ifndef AMBIENT_UTILS_AUXILIARY
#define AMBIENT_UTILS_AUXILIARY

namespace ambient {
    inline void playout()                { ambient::controller.flush();                  }
    inline void conditional_playout()    { ambient::controller.conditional_flush();      }
    inline void mute()                   { ambient::controller.muted = true;             }
    inline void unmute()                 { ambient::controller.muted = false;            }
    inline bool verbose()                { return (rank() ? false : true);               }
}

#endif
