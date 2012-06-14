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
            inline bool requested();
            inline bool trylock();
            inline void unlock();
            inline std::list<cfunctor*>& get_assignments();
            inline std::list<size_t>& get_path();
            void* header;
            void* data;
            bool request;
            bool locked;
            std::list<cfunctor*> assignments;
            std::list<size_t> path;
        };
    
        inline ~layout();
        inline layout(size_t, dim2, dim2);
        inline void embed(void* memory, size_t x, size_t y, size_t bound); // fires controllers::unlock_revision if complete
        inline entry& get(size_t x, size_t y);
        inline void mesh();
        inline size_t id();
        inline const memspec& get_spec() const;
        inline size_t   get_master();
        inline dim2     get_dim() const;
        inline dim2     get_mem_dim() const;
        inline dim2     get_grid_dim() const;
        group* placement;
        size_t master;
        size_t * gid;
        size_t sid;
        dim2   mem_dim;                 // size of distribution blocks
        dim2   mesh_dim;                // size of the grid (reserved) 
        dim2   grid_dim;                // size of the grid 
        dim2   dim;                     // total size of the revision (typed)
        memspec spec;
    private:
        std::vector< layout::entry* > entries;
    }; 
    // }}}

} } }

#endif
