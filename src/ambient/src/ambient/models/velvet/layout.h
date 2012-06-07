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
            inline ~entry();
            inline operator char* (){ return (char*)this->data; }
            inline operator double* (){ return (double*)this->data; }
            inline operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
            inline void set_memory(void* memory, size_t bound);
            inline void* get_memory();
            inline bool valid();
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
        class marker {
        public:
            inline marker();
            inline bool marked(size_t x, size_t y);
            inline void mark(size_t x, size_t y);
            inline void clear();
            bool active;
            size_t xmarker, ymarker;
        };
    
        inline ~layout();
        inline layout(void*, dim2, dim2);
        inline void embed(void* memory, size_t x, size_t y, size_t bound); // fires controllers::unlock_revision if complete
        inline entry* get(size_t x, size_t y);
        inline void mesh();
        inline size_t id();
        inline void*  get_container_type() const;
        inline size_t get_master();
        inline void   set_dim(dim2);
        inline dim2   get_dim() const;
        inline dim2   get_mem_dim() const;
        inline dim2   get_grid_dim() const;
        std::vector< std::vector<layout::entry*> > entries;
        group* placement;
        size_t master;
        size_t * gid;
        size_t sid;
        dim2   mem_dim;                 // size of distribution blocks
        dim2   mesh_dim;                // size of the grid (reserved) 
        dim2   grid_dim;                // size of the grid 
        dim2   dim;                     // total size of the revision (typed)
        void*  container;
    }; 
    // }}}

} } }

#endif
