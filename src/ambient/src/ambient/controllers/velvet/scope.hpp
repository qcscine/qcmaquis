#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { 

    using ambient::controllers::velvet::controller;

    enum scope_t {
        base,
        single,
        shared
    };

    template<scope_t T>
    class scope {};

    template<>
    class scope<base> : public controller::scope {
    public:
        scope(){
            this->round = ambient::channel.volume;
            this->state = ambient::rank() ? REMOTE : LOCAL;
            this->sector = 0;
            /* default common scope:
            this->round = ambient::channel.volume;
            this->state = COMMON;
            this->sector = -1; */
        }
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
            this->state = (this->sector == ambient::rank()) ? LOCAL : REMOTE;
            //printf("switcing scope to %d!\n", this->sector);
        }
       ~scope(){
            if(effect && !--effect) factor = 1;
            ambient::controller.pop_context();
            //printf("reverted to %d!\n", ambient::controller.context->sector);
        }
    };

    template<>
    class scope<shared> : public controller::scope {
    public:
        scope(){
            ambient::controller.uniform = true;
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->state = COMMON;
            this->sector = -1;
            //printf("switcing scope to %d!\n", this->sector);
        }
       ~scope(){
            ambient::controller.uniform = false;
            ambient::controller.pop_context();
            //printf("reverted to %d!\n", ambient::controller.context->sector);
        }
    };
}

#endif
