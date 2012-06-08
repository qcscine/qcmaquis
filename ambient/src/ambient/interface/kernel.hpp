#ifndef AMBIENT_INTERFACE_KERNELS
#define AMBIENT_INTERFACE_KERNELS
#include "ambient/utils/timings.hpp"
 
//#define AMBIENT_PUSH_TIMINGS
#ifdef AMBIENT_PUSH_TIMINGS
    #define __A_TIME(name) static __a_timer time(name); time.begin();
    #define __A_TIME_STOP time.end();
#else
    #define __A_TIME(name) 
    #define __A_TIME_STOP 
#endif

namespace ambient {

    using ambient::controllers::velvet::cfunctor;
    using ambient::models::velvet::revision;

    template<typename FP, FP fp>
    struct kernel_inliner { 
        private:
        static void invoke  (){}
        static void latch   (){}
        static void cleanup (){}
        static void weight  (){}
        static void place   (){}
    };

    #include "ambient/interface/pp/kernel_inliner.pp.hpp"

    template<class K>
    class kernel : public cfunctor
    {
    public:
        virtual ~kernel()         { kernel_inliner<typename K::F,&K::c>::cleanup(this);                            }
        virtual void weight()     { kernel_inliner<typename K::F,&K::c>::weight(this);                             }
        virtual void place()      { kernel_inliner<typename K::F,&K::c>::place(this);                              }
        virtual void computation(){ kernel_inliner<typename K::F,&K::c>::invoke((K*)this); this->check_complete(); }
        virtual void logistics()  { kernel_inliner<typename K::F,&K::l>::invoke((K*)this); this->check_complete(); }
        virtual bool pretend()    { return false; }
        inline void pin(revision& r, int x, int y){ r.content->entries[x][y]->assignments.push_back(this); }
        inline void assign(revision& r, int x, int y){ ambient::controller.ifetch_block(r, x, y); }
        inline void ctxt_select(const char* sql){ this->set_group(channel.world()); }

        inline void block_outright_assign(revision& r){
            size_t sizex = r.content->grid_dim.x;
            size_t sizey = r.content->grid_dim.y;
            for(int x = 0; x < sizex; x++)
                for(int y = 0; y < sizey; y++)
                    this->assign(r, x, y);
        }
        inline void block_outright_pin(revision& r){
            size_t sizex = r.content->grid_dim.x;
            size_t sizey = r.content->grid_dim.y;
            this->add_condition(sizey*sizex);
            for(int x = 0; x < sizex; x++)
                for(int y = 0; y < sizey; y++)
                    this->pin(r, x, y);
            this->block_outright_assign(r);
        }
        inline void block_2d_cycle_assign(revision& r){ 
            this->block_outright_assign(r); 
        }
        inline void block_2d_cycle_pin(revision& r){ 
            this->block_outright_pin(r); 
        }
    };

    template<class K>
    class kernel_atomic : public cfunctor
    {
    public:
        virtual ~kernel_atomic()  { kernel_inliner<typename K::F,&K::c>::cleanup(this);                            }
        virtual void weight()     { kernel_inliner<typename K::F,&K::c>::weight(this);                             }
        virtual void place()      { kernel_inliner<typename K::F,&K::c>::place(this);                              }
        virtual void computation(){ kernel_inliner<typename K::F,&K::c>::invoke((K*)this); this->check_complete(); }
        virtual void logistics()  { kernel_inliner<typename K::F,&K::l>::invoke((K*)this); this->check_complete(); }
        virtual bool pretend()    { return false; }
        inline void ctxt_select(const char* sql){ this->set_group(channel.world()); }
        kernel_atomic(){ 
            this->workload++; 
        }
        inline void pin(revision& r){ 
            r.content->entries[0][0]->assignments.push_back(this);
            ambient::controller.ifetch_block(r,0,0);
        }
        inline void assign(revision& r){ 
            ambient::controller.ifetch_block(r,0,0);
        }
    };

    template<class K>
    class kernel_unpinned : public cfunctor
    {
    public:
        virtual ~kernel_unpinned(){ kernel_inliner<typename K::F,&K::c>::cleanup(this);                            }
        virtual void weight()     { kernel_inliner<typename K::F,&K::c>::weight(this);                             }
        virtual void place()      { kernel_inliner<typename K::F,&K::c>::place(this);                              }
        virtual void computation(){ kernel_inliner<typename K::F,&K::c>::invoke((K*)this); this->check_complete(); }
        inline void pin(revision& r, int x, int y){ r.content->entries[x][y]->assignments.push_back(this); }
        inline void assign(revision& r, int x, int y){ ambient::controller.ifetch_block(r, x, y); }
        inline void ctxt_select(const char* sql){ this->set_group(channel.world()); }
        virtual void logistics(){
            kernel_inliner<typename K::F,&K::l>::invoke((K*)this); 
            this->lock();
            if(this->workload == 1) controller.execute_free_mod(this);
            else this->workload--;
            this->unlock();
        }
        virtual bool pretend(){
            this->lock();
            this->workload--;
            this->unlock();
            if(!this->workload) return false;
            return true;
        }
        inline void atomic_conditional_assign(revision& r){ 
            this->add_condition();
            r.content->entries[0][0]->assignments.push_back(this);
            ambient::controller.ifetch_block(r, 0, 0);
        }
        inline void block_outright_conditional_assign(revision& r){
            size_t sizex = r.content->grid_dim.x;
            size_t sizey = r.content->grid_dim.y;
            this->add_condition(sizey*sizex);
            for(int x = 0; x < sizex; x++)
                for(int y = 0; y < sizey; y++)
                    this->pin(r, x, y);
            for(int x = 0; x < sizex; x++)
                for(int y = 0; y < sizey; y++)
                    this->assign(r, x, y);
        }
        inline void block_2d_cycle_conditional_assign(revision& r){ 
            this->block_outright_conditional_assign(r); 
        }
    };
}

#include "ambient/interface/pp/push.pp.hpp"    
#endif
