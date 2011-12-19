#ifndef AMBIENT_CORE_P_OBJECT_H
#define AMBIENT_CORE_P_OBJECT_H
#include "ambient/core/memblock.h"
#include "ambient/auxiliary.h"
#include <boost/scoped_ptr.hpp>
#include <utility>
#include <list>
#include <vector>

namespace ambient {

    namespace groups { class group; }
    class memblock;
    enum  p_state { ABSTRACT, COMPOSING, GENERIC, PROXY, SERIAL };

    class p_object {
    protected:
        p_object();
    public:
       ~p_object();
// state description knobs /////////////////////
        p_state             state;            //
        bool                consted;          //
////////////////////////////////////////////////

// identification //////////////////////////////
        unsigned int*       group_id;         //
        unsigned int        id;               //
        size_t              timestamp;        //
        groups::group*      scope;            //
        groups::group*      xscope;           //
////////////////////////////////////////////////

// data driven fields //////////////////////////
        size_t              t_size;           //
        p_object*           reference;        // pointer to this profile (this on init - can be changed in proxy objects)
        core::layout_table* layout;           // spatial layout of the profile
        std::vector< std::vector<memblock*> > //
                            skeleton;         //
        size_t              reserved_x;       // skeleton reservation
        size_t              reserved_y;       //
        p_object*           associated_proxy; //
// self-operations /////////////////////////////
        core::operation*    init;             //
        void(*reduce)(memblock*,void*);       //
////////////////////////////////////////////////

// engine fields ///////////////////////////////
        block_packet_t<typename traits::type>*     packet_type;      //
        block_packet_t<typename traits::type>*     xpacket_type;     //
        memblock*           default_block;    //
// lapack knobs ////////////////////////////////
        void*               data;             // 
        size_t              solid_lda;        //
////////////////////////////////////////////////

// spatial configuration ///////////////////////
        dim2                dim;              //
    private:                                  //
        dim2                work_dim;         // work-item size of distribution blocks
        dim2                mem_dim;          // work-item size of cpu streaming multiprocessor workload fractions
        dim2                item_dim;         // size of work-item (i.e. 128) 
        dim2                gpu_dim;          // work-item size of gpgpu smp workload fractions
////////////////////////////////////////////////
    public:
        p_object&          operator >>(dim2 mem_dim);
        p_object&          operator , (dim2 dim);

        void                constant();
        void                inconstant();
        bool                is_associated_proxy();

        void                set_init(core::operation* op);
        core::operation*    get_init() const;
        void                associate_proxy(void(*R)(memblock*,void*));
        void                reblock();
        void                solidify(std::vector<core::layout_table::entry> entries);
        void                disperse(std::vector<core::layout_table::entry> entries);

        std::pair<unsigned int*,size_t> 
                            get_id();
        groups::group*      get_xscope();
        groups::group*      get_scope();
        void                set_scope(groups::group* scope);
        bool                involved();
        bool                xinvolved();

        void                touch();
        void                preprocess();
        void                postprocess();             // resets init marker
        void                postprocess(int i, int j); // proceed with necessary memory allocations
        void                finalize();                // proceed with proxy updates (various reduces)
        void                clean();                   // cleanup for proxy/layout junk

        dim2                get_block_id();
        void                set_default_block(int i, int j = 0);
        size_t              get_block_lda();

        memblock*           block(int i, int j = 0) const;
        memblock&           operator()(int i, int j = 0);
        void*               get_data();
        template<typename T> operator T*()
        { return (T*)this->get_data();   }
    private:
        template<typename T> operator T ();
    public:
// parameters can be set specifically for the profile
        dim2                get_dim()       const;
        void                set_dim(dim2 dim);
        dim2                get_work_dim()  const;
        void                set_work_dim(dim2 dim);
        dim2                get_gpu_dim()   const;
        void                set_gpu_dim(dim2 dim);
        dim2                get_grid_dim()  const;
        dim2                get_mem_t_dim() const;
        dim2                get_mem_dim()   const;
        void                set_mem_dim(dim2 dim);
        dim2                get_item_dim()  const;
        void                set_item_dim(dim2 dim);
        size_t              get_bound()     const;
    };

    p_object& operator>>(p_object* instance, dim2 mem_dim);
    void accept_block(groups::packet_manager::typed_q& in_q);
}
#endif
