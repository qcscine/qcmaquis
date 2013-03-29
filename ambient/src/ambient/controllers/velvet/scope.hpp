#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { 

    using ambient::controllers::velvet::controller;

    template<scope_t T>
    class scope {};

    template<>
    class scope<base> : public controller::tunable_scope {
    public:
        scope(){
            this->tunable = true;
            this->round = ambient::channel.volume;
            this->state = ambient::rank() ? ambient::remote : ambient::local;
            this->sector = 0;
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
            this->tunable = false;
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->sector = (++iterator %= round*factor)/factor;
            this->state = (this->sector == ambient::rank()) ? ambient::local : ambient::remote;
        }
       ~scope(){
            if(effect && !--effect) factor = 1;
            ambient::controller.pop_context();
        }
    };

    template<>
    class scope<shared> : public controller::scope {
    public:
        scope(){
            this->tunable = false;
            ambient::controller.set_context(this);
            this->round = ambient::channel.volume;
            this->state = ambient::common;
        }
       ~scope(){
            ambient::controller.pop_context();
        }
    };
}

#endif
