#ifndef AMBIENT_CORE_P_PROFILE_H
#define AMBIENT_CORE_P_PROFILE_H
#include "ambient/core/memblock.h"
#include "ambient/auxiliary.h"
#include <boost/scoped_ptr.hpp>
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    namespace groups { class group; }
    class memblock;
    enum  p_state { ABSTRACT, COMPOSING, GENERIC, PROXY };

    class p_profile {
    protected:
        p_profile();
    public:
        p_state             state;
        unsigned int*       group_id;
        unsigned int        id;
        size_t              timestamp;
        bool                consted;
        bool                finalized;
        std::pair<int,int>  master_relay;
        p_profile*          profile; // pointer to this profile (this on init - can be changed in proxy objects)
        void*               data;  // pointer to the actual data
        size_t              lda;  // process individual lda
        size_t              solid_lda;  // process solid state matrix lda
        size_t              block_lda;
        groups::group*      scope;
        groups::group*      xscope;
        core::layout_table* layout;
        size_t              reserved_x;
        size_t              reserved_y;
        size_t              t_size;
        dim2                dim;
        block_packet_t*     packet_type;
        block_packet_t*     xpacket_type;
        std::vector< std::vector<memblock*> > skeleton;
        memblock*          default_block;
        void(*init)(memblock*);
        void(*reduce)(memblock*,void*);
        p_profile*          associated_proxy;
    private:
        bool                valid;
        dim2                distr_dim;   // work-item size of distribution blocks
        dim2                mem_dim;   // work-item size of cpu streaming multiprocessor workload fractions
        dim2                item_dim;    // size of work-item (i.e. 128) 
        dim2                gpu_dim;     // work-item size of gpgpu smp workload fractions
    public:
        void operator=(const p_profile& profile);
        p_profile & operator>>(dim2 dim_distr);
        p_profile & operator,(dim2 dim);

        void constant();
        void inconstant();
        p_profile* associate_proxy(p_profile* proxy, void(*R)(memblock*,void*));

        bool is_proxy();
        void reblock();
        void set_id(std::pair<unsigned int*,size_t> group_id);
        std::pair<unsigned int*,size_t> get_id();
        void set_master(int master);
        int get_master();
        int get_xmaster();
        void set_default_block(int i, int j = 0);
        dim2 get_block_id();

        p_profile* dereference(); // finds out if the profile pointer is up to date
        void touch();
        void preprocess();
        void postprocess();             // resets init marker
        void postprocess(int i, int j); // proceed with necessary memory allocations
        void finalize();    // proceed with proxy updates (various reduces)
        void clean();       // cleanup for proxy/layout junk
        size_t get_block_lda();

        void set_scope(groups::group* scope);
        groups::group* get_scope();
        groups::group* get_xscope();
        bool involved();
        bool xinvolved();

        memblock* block(int i, int j = 0) const;
        memblock& operator()(int i, int j = 0);

// parameters can be set specifically for the profile
        dim2 get_dim()         const;
        dim2 get_distr_dim()   const;
        dim2 get_gpu_dim()     const;
        dim2 get_grid_dim()    const;
        dim2 get_mem_dim()   const;
        dim2 get_mem_t_dim() const;
        dim2 get_item_dim()    const;
        void(*get_init() const)(memblock*);

        void imitate(p_profile* profile);
        void solidify(std::vector<core::layout_table::entry> entries);
        void disperse(std::vector<core::layout_table::entry> entries);

    private:
        template<typename T> operator T ();
    public:
        template<typename T> operator T* ()
        { return (T*)this->get_data();    }
        size_t get_bound() const;
        void* get_data();
        void set_dim(dim2 dim);
        void set_distr_dim(dim2 dim);
        void set_gpu_dim(dim2 dim);
        void set_mem_dim(dim2 dim);
        void set_item_dim(dim2 dim);
        void set_init(void(*)(memblock*));
        void invalidate();
        bool is_valid();
    };

    p_profile& operator>>(p_profile* instance, dim2 dim_distr);
    void accept_block(groups::packet_manager::typed_q& in_q);
}
#endif
