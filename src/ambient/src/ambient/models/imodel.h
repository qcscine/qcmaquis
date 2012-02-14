#ifndef AMBIENT_INTERFACE_MODEL_H
#define AMBIENT_INTERFACE_MODEL_H
#include "ambient/ambient.h"
#include "ambient/utils/dim2.h"
//#include "ambient/utils/delegate.hpp"

namespace ambient { namespace models {

    class imodel {
    public:
        class object;
        class modifier {
            public:
            virtual void invoke() = 0;
            virtual void weight() = 0;
            virtual size_t get_weight() = 0;
            virtual void set_weight(size_t) = 0;
            virtual void set_group(channels::group* grp) = 0;
            virtual void set_vellum(object&) = 0;
            virtual object& get_vellum() = 0;
        };
        class revision;
        class layout {
            public:
            class entry {
            public:
                virtual bool valid() = 0;
                virtual operator double* () = 0;
                virtual operator std::complex<double>* () = 0;
                virtual void* get_memory() = 0;
                virtual void set_memory(void* memory, size_t bound) = 0;
            };
            class path {
            public:
                virtual void push_node(size_t n) = 0;
                virtual size_t pop_node() = 0;
                virtual bool empty() = 0;
            };
            virtual void push_path(size_t i, size_t j, size_t node) = 0;
            virtual size_t pop_path(size_t i, size_t j) = 0;
            virtual path* get_path(size_t i, size_t j) = 0;
            virtual entry* get(size_t i, size_t j) = 0;
            virtual void set_dim(dim2) = 0;
            virtual dim2 get_dim() const = 0;
            virtual dim2 get_grid_dim() const = 0;
            virtual dim2 get_mem_grid_dim() const = 0;
            virtual dim2 get_mem_dim() const = 0;
            virtual dim2 get_work_dim() const = 0;
            virtual dim2 get_item_dim() const = 0;
            virtual void set_revision(revision* r) = 0;
            virtual std::pair<size_t*,size_t> id() = 0;
            virtual size_t get_mem_size() const = 0;
            virtual size_t get_mem_lda() const = 0;
            virtual size_t get_master() = 0;
        };
        class revision {
        public:
            virtual void set_dim(dim2) = 0;
            virtual std::pair<size_t*,size_t> id() = 0;
            virtual imodel::layout::entry& operator()(size_t i, size_t j) = 0;
            virtual void add_modifier(modifier* m) = 0;
            virtual imodel::layout& get_layout() = 0;
            virtual channels::group* get_placement() = 0;
            virtual void reduce(void(*)(void*,void*)) = 0;
        };
        class object {
        public:
            virtual void add_revision(imodel::layout* l) = 0;
            virtual imodel::revision& revision(size_t offset) const = 0;
            virtual dim2 get_dim() const = 0;
            virtual size_t get_t_size() const = 0;
            virtual channels::group* get_placement() = 0;
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
