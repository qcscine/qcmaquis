#ifndef AMBIENT_MODELS_V_MODEL_H
#define AMBIENT_MODELS_V_MODEL_H
#include "ambient/ambient.h"
#include "ambient/models/v_model.h"
#include "ambient/utils/singleton.hpp"
#include "ambient/utils/hashmap.h"
#include "ambient/utils/dim2.h"
#include <list>

extern pthread_key_t pthread_tid;

#define GET_TID 0 //*(size_t*)pthread_getspecific(pthread_tid)

namespace ambient { namespace models {
// fires:
// object revision ready
// (to controller)

    class v_model : public singleton< v_model > {
    public:
        // {{{ layout memory-specific model
        class revision;
        class object;
        class reduction;
        class modifier;
        class layout {
        public:
            class entry {
            public:
                entry();
               ~entry();
                entry(void*, size_t);
                inline operator char* (){ return (char*)this->data; }
                inline operator double* (){ return (double*)this->data; }
                inline operator std::complex<double>* (){ return (std::complex<double>*)this->data; }
                void set_memory(void* memory, size_t bound);
                void* get_memory();
                bool valid();
                bool requested();
                bool trylock();
                void unlock();
                std::list<models::v_model::modifier*>& get_assignments();
                std::list<size_t>& get_path();
                void* header;
                void* data;
                bool request;
                bool locked;
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
            void set_revision(v_model::revision* r);
            std::pair<size_t*,size_t> id();
            size_t get_master();
            size_t get_mem_size() const;
            size_t get_mem_lda() const;

            layout& operator>>(dim2);       // set mem_dim
            layout& operator, (dim2);       // set item_dim
       
            void   set_dimensions(dim2,dim2);
            void   set_dim(dim2);
            dim2   get_dim() const;
            dim2   get_mem_dim() const;
            dim2   get_item_dim() const;
            dim2   get_grid_dim() const;
            std::vector< std::vector<v_model::layout::entry*> > entries;
            v_model::layout::marker marker;
            v_model::revision* revision;
            size_t master;
            size_t * gid;
            size_t sid;
            size_t t_size;
            dim2   mem_dim;                 // size of distribution blocks
            dim2   item_dim;                // size of work-item inside workgroup
            dim2   mesh_dim;                // size of the grid (reserved) 
            dim2   grid_dim;                // size of the grid 
            dim2   dim;                     // total size of the revision (typed)
        }; 
        // }}}
        // {{{ revision property (modifier)
        class modifier {
        protected:
            public:
            virtual ~modifier(){};
            virtual void invoke() = 0;
            virtual void weight() = 0;
            virtual size_t get_weight() = 0;
            virtual void set_weight(size_t) = 0;
            virtual void set_group(channels::group* grp) = 0;
            virtual channels::group* get_group() = 0;
            virtual void set_vellum(revision&) = 0;
            virtual revision& get_vellum() = 0;
            virtual revision* get_pin() = 0;
            virtual bool pretend() = 0;
            virtual void add_condition() = 0;
        };
        // }}}
        // {{{ revisioned object impl-specific model
        class revision
        {
        public:
            inline v_model::layout::entry& operator()(size_t i, size_t j){
                v_model::layout::entry* e = this->layout->entries[i][j];
                if(e->header == NULL) return this->alloc_block(i, j); // serial
                return *e;
            } 
            inline v_model::layout& get_layout(){
                return *this->layout;
            }
            v_model::layout::entry& alloc_block(size_t i, size_t j);
            revision(v_model::object*, v_model::layout*);
           ~revision();
            v_model::layout::entry* block(size_t i, size_t j = 0);
            void add_modifier(v_model::modifier* m);
            std::list<v_model::modifier*>& get_modifiers();
            std::pair<size_t*,size_t> id();
            v_model::object& get_object();
            channels::group* get_placement();
            void set_placement(channels::group*);
            v_model::reduction* get_reduction();
            void set_reduction();
            v_model::modifier* get_generator();
            void set_generator(v_model::modifier*);
            void set_dim(dim2);
            dim2 get_dim();
            size_t number;
            v_model::object* const object;
            v_model::layout* const layout;
            channels::group* placement;
            v_model::modifier* generator;
            v_model::reduction* reduction;
            std::list<v_model::modifier*> modifiers;
        };
        class fast_revision // evil... that's evil
        {
        public:
            inline v_model::layout::entry& operator()(size_t i, size_t j){
                // while(!this->valid()) pthread_yield(); // can be done if MPI enabled
                return *((revision*)this)->layout->entries[i][j];
            }
            inline v_model::layout& get_layout(){
                return *((revision*)this)->layout;
            }
        };
        class reduction
        {
        public:
            class reductionq {
            public:
                reductionq();
                void push(v_model::layout::entry*);
            };
            reduction(v_model::revision*);
           ~reduction();
            v_model::layout::entry* block(size_t i, size_t j);
            v_model::layout::entry& operator()(size_t i, size_t j);
            std::vector< std::vector<reductionq*> > entries;
            v_model::revision* revision;
        };
        class object 
        {            // revision tracking mechanism (target selector)
        protected:
            object();
        public:
           ~object();
            void add_revision(v_model::layout* l);
            v_model::revision& revision(size_t offset) const;

            // light-weight user-interface //
            inline v_model::fast_revision& ui_c_revision_0() const { 
                return *(fast_revision*)this->revisions[this->thread_revision_base[GET_TID]]; 
            }
            inline v_model::revision& ui_c_revision_1() const { 
                return *this->revisions[this->thread_revision_base[GET_TID] + 1]; 
            }
            // light-weight user-interface //

            dim2 get_dim() const;
            void set_dim(dim2);
            size_t get_t_size() const;
            size_t get_revision_base() const;
            size_t get_thread_revision_base() const;
            void set_revision_base(size_t);
            void set_thread_revision_base(size_t);
            std::vector<v_model::revision*> revisions;
            size_t t_size;
            size_t* thread_revision_base;
            size_t revision_base;
            dim2 init_dim;
        };
        // }}}
    public: 
        v_model();
        void add_revision(v_model::object* obj);
        void update_revision(v_model::revision* r, channels::group* placement);
        v_model::revision* get_revision(size_t* hash, size_t hash_len, size_t id) const;
        v_model& operator>>(dim2);
        v_model& operator, (dim2);
       ~v_model();
    private:
        hashmap map;
        dim2 mem_dim;
        dim2 item_dim;
    };

    // free functions for mangling the data {{{ 
    void* solidify(const v_model::object& o);
    void  disperse(void* data, v_model::object& o);
    // }}}
} }

namespace ambient {
    extern models::v_model& model;
}

#endif
