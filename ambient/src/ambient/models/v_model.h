#ifndef AMBIENT_VELVET_MODEL_H
#define AMBIENT_VELVET_MODEL_H
#include "ambient/ambient.h"
#include "ambient/models/imodel.h"
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/hashmap.h"
#include "ambient/utils/dim2.h"
#include <list>

namespace ambient { namespace models {
// fires:
// object revision ready
// (to controller)

    class v_model : public imodel, public singleton< v_model > {
    public:
        // {{{ layout memory-specific model
        class revision;
        class layout : public imodel::layout {
        public:
            class entry : public imodel::layout::entry {
            public:
                entry();
                operator double* ();
                operator std::complex<double>* ();
                void set_memory(void* memory, size_t bound);
                void* get_memory();
                bool valid();
                void* header;
                void* data;
            };
            class path : public imodel::layout::path {
            public:
                path();
                void add_node(size_t n);
                std::list<size_t> route;
            };
            class marker {
            public:
                marker();
                bool marked(size_t i, size_t j);
                void mark(size_t i, size_t j);
                void clear();
                bool active;
                size_t imarker, jmarker;
            };
        
           ~layout();
            layout(dim2,size_t);
            void embed(void* memory, size_t i, size_t j); // fires controllers::unlock_revision if complete
            void add_path(size_t i, size_t j, size_t node);
            entry* get(size_t i, size_t j);
            void mesh();
            void set_revision(imodel::revision* r);
            std::pair<size_t*,size_t> id();
            size_t get_master();
            size_t get_mem_size() const;

            layout& operator>>(dim2);       // set mem_dim
            layout& operator, (dim2);       // set work_dim
        
            void   add_route(size_t node);
            void   set_dim(dim2);
            dim2   get_dim() const;
            dim2   get_mem_dim() const;
            dim2   get_work_dim() const;
            dim2   get_item_dim() const;
            dim2   get_grid_dim() const;
            dim2   get_mem_grid_dim() const;
            dim2   mem_dim;                 // size of distribution blocks
            dim2   work_dim;                // size of csmp workload groups
            dim2   item_dim;                // size of work-item inside workgroup
            dim2   mesh_dim;                // size of the grid (reserved) 
            dim2   dim;                     // total size of the revision (typed)
            std::vector< std::vector<v_model::layout::entry*> > entries;
            std::vector< std::vector<v_model::layout::path*> > paths;
            v_model::layout::marker marker;
            v_model::revision* revision;
            size_t master;
            size_t * gid;
            size_t sid;
            size_t t_size;
        }; 
        // }}}
        // {{{ revision property (modifier)
        class modifier : public imodel::modifier {
        protected:
            modifier();
            size_t locks;
            size_t workload;
        };
        // }}}
        // {{{ revisioned object impl-specific model
        class revision : public imodel::revision
        {
        public:
            revision(imodel::layout*);
           ~revision();
            imodel::layout::entry* block(size_t i, size_t j = 0);
            imodel::layout::entry& operator()(size_t i, size_t j);
            void add_modifier(imodel::modifier* m);
            std::pair<size_t*,size_t> id();
            imodel::layout& get_layout();
            channels::group* get_placement();
            void   set_dim(dim2);
            size_t number;
            imodel::layout* const layout;
            channels::group* placement;
            std::list<imodel::modifier*> modifiers;
        };
        class object: public imodel::object
        {            // revision tracking mechanism (target selector)
        protected:
            object();
        public:
           ~object();
            void add_revision(imodel::layout* l);
            v_model::revision& revision(size_t offset) const;
            dim2 get_dim() const;
            size_t get_t_size() const;
            channels::group* get_placement();
            std::vector<v_model::revision*> revisions;
            size_t t_size;
            size_t base; // revision base ??
            channels::group* placement;
            dim2 dim;
        };
        // }}}
    public: 
        v_model();
        void add_revision(imodel::object* obj);
        void update_revision(channels::group* placement, imodel::revision* r);
        v_model::revision* get_revision(size_t* hash, size_t hash_len, size_t id) const;
        v_model& operator>>(dim2);
        v_model& operator, (dim2);
        dim2 get_mem_dim() const;
        dim2 get_work_dim() const;
        dim2 get_item_dim() const;
       ~v_model();
    private:
        hashmap map;
        dim2 mem_dim;
        dim2 work_dim;
        dim2 item_dim;
    };

    // free functions for mangling the data {{{ 
    void* solidify(imodel::revision& instance);
    void  disperse(void* data, imodel::revision& instance);
    void  reduce(v_model::modifier* r);
    void  init(v_model::modifier* i);
    // }}}
} }
#endif
