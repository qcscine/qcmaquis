#ifndef AMBIENT_INTERFACE_MODEL
#define AMBIENT_INTERFACE_MODEL
#include "ambient/ambient.h"
#include "ambient/utils/dim2.h"
//#include "ambient/utils/delegate.hpp"
#include <list>

namespace ambient { namespace models {

    class imodel {
    public:
        class object;
        class revision;
        class reduction;
        class modifier {
            public:
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
        class layout {
            public:
            class entry {
            public:
                virtual bool trylock() = 0;
                virtual void unlock() = 0;
                virtual bool valid() = 0;
                virtual bool requested() = 0;
                virtual operator char* () = 0;
                virtual operator double* () = 0;
                virtual operator std::complex<double>* () = 0;
                virtual void* get_memory() = 0;
                virtual void set_memory(void* memory, size_t bound) = 0;
                virtual std::list<modifier*>& get_assignments() = 0;
                virtual std::list<size_t>& get_path() = 0;
            };
            virtual void mark(size_t i, size_t j) = 0;
            virtual bool marked(size_t i, size_t j) = 0;
            virtual void embed(void* memory, size_t i, size_t j, size_t bound) = 0;
            virtual entry* get(size_t i, size_t j) = 0;
            virtual void set_dimensions(dim2,dim2) = 0;
            virtual void set_dim(dim2) = 0;
            virtual dim2 get_dim() const = 0;
            virtual dim2 get_grid_dim() const = 0;
            virtual dim2 get_mem_grid_dim() const = 0;
            virtual dim2 get_mem_dim() const = 0;
            virtual dim2 get_item_dim() const = 0;
            virtual void set_revision(revision* r) = 0;
            virtual std::pair<size_t*,size_t> id() = 0;
            virtual size_t get_mem_size() const = 0;
            virtual size_t get_mem_lda() const = 0;
            virtual size_t get_master() = 0;
        };
        class revision {
        public:
            typedef void(*voidfp)();
            virtual void set_dim(dim2) = 0;
            virtual dim2 get_dim() = 0;
            virtual std::pair<size_t*,size_t> id() = 0;
            virtual imodel::layout::entry& operator()(size_t i, size_t j) = 0;
            virtual imodel::layout::entry* block(size_t i, size_t j) = 0;
            virtual void add_modifier(modifier* m) = 0;
            virtual std::list<modifier*>& get_modifiers() = 0;
            virtual imodel::object& get_object() = 0;
            virtual imodel::layout& get_layout() = 0;
            virtual channels::group* get_placement() = 0;
            virtual void set_placement(channels::group*) = 0;
            virtual modifier* get_generator() = 0;
            virtual void set_generator(modifier*) = 0;
            virtual reduction* get_reduction() = 0;
            virtual void set_reduction() = 0;
        };
        class reduction {
        public:
            virtual imodel::layout::entry& operator()(size_t i, size_t j) = 0;
        };
        class object {
        public:
            virtual void add_revision(imodel::layout* l) = 0;
            virtual imodel::revision& revision(size_t offset) const = 0;
            virtual dim2 get_dim() const = 0;
            virtual void set_dim(dim2) = 0;
            virtual size_t get_t_size() const = 0;
            virtual size_t get_revision_base() const = 0;
            virtual size_t get_thread_revision_base() const = 0;
            virtual void set_revision_base(size_t) = 0;
            virtual void set_thread_revision_base(size_t) = 0;
        };
        virtual void add_revision(object* obj) = 0;
        virtual void update_revision(revision* r, channels::group* placement) = 0;
        virtual revision* get_revision(size_t* hash, size_t hash_len, size_t id) const = 0;
        virtual imodel& operator>>(dim2) = 0;
        virtual imodel& operator, (dim2) = 0;
        //delegate object_modified;
        //delegate object_released;
    };


} }

namespace ambient {
    extern models::imodel&  model;
}

#endif
