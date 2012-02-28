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
                bool requested();
                std::list<models::imodel::modifier*>& get_assignments();
                std::list<size_t>& get_path();
                void* header;
                void* data;
                bool request;
                std::list<modifier*> assignments;
                std::list<size_t> path;
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
            void mark(size_t i, size_t j);
            bool marked(size_t i, size_t j);
            void embed(void* memory, size_t i, size_t j, size_t bound); // fires controllers::unlock_revision if complete
            entry* get(size_t i, size_t j);
            void mesh();
            void set_revision(imodel::revision* r);
            std::pair<size_t*,size_t> id();
            size_t get_master();
            size_t get_mem_size() const;
            size_t get_mem_lda() const;

            layout& operator>>(dim2);       // set mem_dim
            layout& operator, (dim2);       // set work_dim
        
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
            typedef void(*voidfp)();
            revision(imodel::object*, imodel::layout*);
           ~revision();
            imodel::layout::entry* block(size_t i, size_t j = 0);
            imodel::layout::entry& operator()(size_t i, size_t j);
            void add_modifier(imodel::modifier* m);
            std::list<imodel::modifier*>& get_modifiers();
            std::pair<size_t*,size_t> id();
            imodel::object& get_object();
            imodel::layout& get_layout();
            channels::group* get_placement();
            void set_placement(channels::group*);
            imodel::modifier* get_generator();
            void set_generator(imodel::modifier*);
            void reduce(voidfp);
            void init(voidfp);
            voidfp get_reduce();
            voidfp get_init();
            void set_dim(dim2);
            size_t number;
            imodel::object* const object;
            imodel::layout* const layout;
            channels::group* placement;
            imodel::modifier* generator;
            std::list<imodel::modifier*> modifiers;
            voidfp reduction;
            voidfp initialization;
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
            size_t get_revision_base() const;
            size_t get_thread_revision_base() const;
            void set_revision_base(size_t);
            void set_thread_revision_base(size_t);
            std::vector<v_model::revision*> revisions;
            size_t t_size;
            pthread_key_t thread_revision_base;
            size_t revision_base;
            dim2 dim;
        };
        // }}}
    public: 
        v_model();
        void add_revision(imodel::object* obj);
        void update_revision(imodel::revision* r, channels::group* placement);
        v_model::revision* get_revision(size_t* hash, size_t hash_len, size_t id) const;
        v_model& operator>>(dim2);
        v_model& operator, (dim2);
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
    // }}}
} }
#endif
