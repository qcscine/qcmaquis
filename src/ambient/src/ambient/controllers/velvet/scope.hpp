#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { 

    using ambient::controllers::velvet::controller;

    template<scope_t T>
    class scope {};

    template<>
    class scope<base> : public controller::scope {
    public:
        scope(){
            this->round = ambient::channel.volume;
            this->state = ambient::rank() ? ambient::remote : ambient::local;
            this->sector = 0;
            this->factor = AMBIENT_SCOPE_SWITCH_FACTOR;
            this->op_alloc = 0;
            this->op_transfer = 0;
        }
        virtual bool tunable(){ 
            return false; // can be enabled but the algorithm should be finalized
        }
        virtual void consider_allocation(size_t size){
            this->op_alloc += size;
        }
        virtual void consider_transfer(size_t size, ambient::locality l){
            if(l == ambient::common) return;
            this->op_transfer += size;
        }
        virtual void toss(){
            if(this->op_transfer < this->op_alloc){
                if(this->factor < this->op_alloc){
                    this->factor = AMBIENT_SCOPE_SWITCH_FACTOR;
                    ++this->sector %= this->round;
                    this->state = (this->sector == ambient::rank()) ? 
                                  ambient::local : ambient::remote;
                }else{
                    this->factor -= this->op_alloc;
                }
            }
            this->op_alloc = 0;
            this->op_transfer = 0;
        }
        size_t factor;
        size_t op_alloc;
        size_t op_transfer;
    };

    template<>
    class scope<single> : public controller::scope {
    public:
        static int iterator;
        static int factor;
        static int effect;
        static void compact(size_t n){
            iterator = 0;
            if(n <= ambient::channel.volume) return;
            factor = (int)(n / ambient::channel.volume); // iterations before switch
            effect = (int)n;
        }
        scope(){
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->sector = (++iterator %= round*factor)/factor;
            this->state = (this->sector == ambient::rank()) ? ambient::local : ambient::remote;
        }
       ~scope(){
            if(effect && !--effect) factor = 1;
            ambient::controller.pop_context();
        }
        virtual bool tunable(){ 
            return false; 
        }
    };

    template<>
    class scope<shared> : public controller::scope {
    public:
        scope(){
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->state = ambient::common;
        }
       ~scope(){
            ambient::controller.pop_context();
        }
        virtual bool tunable(){ 
            return false; 
        }
    };
}

#endif
