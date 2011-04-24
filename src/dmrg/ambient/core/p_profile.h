#ifndef AMBIENT_CORE_P_PROFILE_H
#define AMBIENT_CORE_P_PROFILE_H
#include "ambient/core/workgroup.h"
#include "ambient/auxiliary.h"
#include <boost/scoped_ptr.hpp>
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    namespace groups { class group; }
    class workgroup;
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
        size_t              group_lda;
        groups::group*      scope;
        groups::group*      xscope;
        core::layout_table* layout;
        size_t              reserved_x;
        size_t              reserved_y;
        size_t              t_size;
        dim3                dim;
        block_packet_t*     packet_type;
        block_packet_t*     xpacket_type;
        std::vector< std::vector<workgroup*> > skeleton;
        workgroup*          default_group;
        void(*init)(workgroup*);
        void(*reduce)(workgroup*,void*);
        p_profile*          associated_proxy;
    private:
        bool                valid;
        dim3                distr_dim;   // work-item size of distribution blocks
        dim3                group_dim;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3                item_dim;    // size of work-item (i.e. 128) 
        dim3                gpu_dim;     // work-item size of gpgpu smp workload fractions
    public:
        void operator=(const p_profile& profile);
        p_profile & operator>>(dim3 dim_distr);
        p_profile & operator,(dim3 dim);

        void constant();
        void inconstant();
        p_profile* associate_proxy(p_profile* proxy, void(*R)(workgroup*,void*));

        bool is_proxy();
        void regroup();
        void set_id(std::pair<unsigned int*,size_t> group_id);
        std::pair<unsigned int*,size_t> get_id();
        void set_master(int master);
        int get_master();
        int get_xmaster();
        void set_default_group(int i, int j = 0, int k = 0);
        dim3 get_group_id();

        p_profile* dereference(); // finds out if the profile pointer is up to date
        void touch();
        void preprocess();
        void postprocess();             // resets init marker
        void postprocess(int i, int j); // proceed with necessary memory allocations
        void finalize();    // proceed with proxy updates (various reduces)
        void clean();       // cleanup for proxy/layout junk
        size_t get_group_lda();

        void set_scope(groups::group* scope);
        groups::group* get_scope();
        groups::group* get_xscope();
        bool involved();
        bool xinvolved();

        workgroup* group(int i, int j = 0, int k = 0) const;
        workgroup& operator()(int i, int j = 0, int k = 0);

// parameters can be set specifically for the profile
        dim3 get_dim()         const;
        dim3 get_distr_dim()   const;
        dim3 get_gpu_dim()     const;
        dim3 get_grid_dim()    const;
        dim3 get_group_dim()   const;
        dim3 get_group_t_dim() const;
        dim3 get_item_dim()    const;
        void(*get_init() const)(workgroup*);

        void imitate(p_profile* profile);
        void solidify(std::vector<core::layout_table::entry>& entries);
        void disperse(std::vector<core::layout_table::entry>& entries);

    private:
        template<typename T> operator T ();
    public:
        template<typename T> operator T* ()
        { return (T*)this->get_data();    }
        size_t get_bound() const;
        void* get_data();
        void set_dim(dim3 dim);
        void set_distr_dim(dim3 dim);
        void set_gpu_dim(dim3 dim);
        void set_group_dim(dim3 dim);
        void set_item_dim(dim3 dim);
        void set_init(void(*)(workgroup*));
        void invalidate();
        bool is_valid();
    };

    p_profile& operator>>(p_profile* instance, dim3 dim_distr);
    void accept_block(groups::packet_manager::typed_q& in_q);
}
#endif
