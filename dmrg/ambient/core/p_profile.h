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

    class p_profile {
    protected:
        p_profile();
    public:
        unsigned int* group_id;
        unsigned int id;
        p_profile & operator>>(dim3 dim_distr);
        p_profile & operator,(dim3 dim);

        void constant();
        void inconstant();
        bool consted;

        void regroup();
        void set_id(std::pair<unsigned int*,size_t> group_id);
        void set_master(int master);
        int get_master();
        int get_xmaster();
        void set_default_group(int i, int j = 0, int k = 0);
        dim3 get_group_id();

        std::pair<int,int> master_relay;
        p_profile* profile; // pointer to this profile (this on init - can be changed in proxy objects)
        p_profile* dereference(); // finds out if the profile pointer is up to date
        void touch();
        void preprocess();
        void postprocess(); // proceed with necessary memory allocations
        void* framework; // includes ambient boundle
        void* data;  // pointer to the actual data
        size_t lda;  // process individual lda
        size_t solid_lda;  // process solid state matrix lda
        size_t group_lda;
        size_t get_group_lda();

        groups::group* scope;
        groups::group* xscope;

        void set_scope(groups::group* scope);
        groups::group* get_scope();
        groups::group* get_xscope();
        bool involved();
        bool xinvolved();

        core::layout_table* layout;

        size_t reserved_x;
        size_t reserved_y;

        const char* type;
        size_t type_size;
        bool proxy;
        
        dim3 dim;
        int** owners;

        std::vector< std::vector<workgroup*> > skeleton;
        workgroup* default_group;

        void(*init_fp)(workgroup* grp);
        workgroup* group(int i, int j = 0, int k = 0) const;
        workgroup& operator()(int i, int j = 0, int k = 0) const;
        block_packet_t* packet_type;

// parameters can be set specifically for the profile
        dim3 get_dim()       const;
        dim3 get_distr_dim() const;
        dim3 get_gpu_dim()   const;
        dim3 get_grid_dim()  const;
        dim3 get_group_dim() const;
        dim3 get_item_dim()  const;

        void imitate(p_profile* profile);
        void solidify();
        void disperse();

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
        void invalidate();
        bool is_valid();
        bool is_inited();
        bool need_init;
        bool preprocessed;
    private:
        bool inited;
        bool valid;
        dim3 distr_dim;   // work-item size of distribution blocks
        dim3 group_dim;   // work-item size of cpu streaming multiprocessor workload fractions
        dim3 item_dim;    // size of work-item (i.e. 128) 
        dim3 gpu_dim;     // work-item size of gpgpu smp workload fractions
    };

    p_profile& operator>>(p_profile* instance, dim3 dim_distr);
    void integrate_block(groups::packet_manager::typed_q& in_q);
}
#endif
