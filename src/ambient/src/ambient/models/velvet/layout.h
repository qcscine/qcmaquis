#ifndef AMBIENT_MODELS_VELVET_LAYOUT
#define AMBIENT_MODELS_VELVET_LAYOUT
#include "ambient/channels/mpi/groups/group.h"

namespace ambient { namespace controllers { namespace velvet {

    class cfunctor;

} } }

namespace ambient { namespace models { namespace velvet {

    using ambient::controllers::velvet::cfunctor;
    using ambient::channels::mpi::group;

    // {{{ layout memory-specific model
    class revision;
    class layout {
    public:
        class entry {
        public:
            inline entry();
            inline operator char* (){ return (char*)this->data; }
            inline operator double* (){ return (double*)this->data; }
            inline operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
            inline void swap(entry&);
            inline void set_memory(void* memory, size_t bound);
            inline void* get_memory();
            inline bool valid();
            inline bool occupied();
            //inline bool requested();
            //inline bool trylock();
            //inline void unlock();
            inline std::list<cfunctor*>& get_assignments();
            //inline std::list<size_t>& get_path();
            void* header;
            void* data;
            //bool request;
            //bool locked;
            std::list<cfunctor*> assignments;
            //std::list<size_t> path;
        };
    
        inline ~layout();
        inline layout(memspec* spec, bool clean = false);
        inline void embed(void* memory, size_t x, size_t y, size_t bound); // fires controllers::unlock_revision if complete
        inline entry& get(size_t x, size_t y);
        //inline size_t get_master();
        //group* placement;
        //size_t master;
        size_t sid;
        size_t lda;
        memspec* spec;
        bool clean;
    private:
        std::vector< layout::entry* > entries;
    }; 
    // }}}

} } }

#endif
