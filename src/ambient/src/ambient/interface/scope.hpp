#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { 

    enum scope_t {
        single,
        shared,
        common
    };

    template<scope_t T>
    class scope {
    public:
    /*    scope(){
            printf("Created Ambient shared scope!\n");
        }
        scope(size_t req){
            printf("Created Ambient shared scope of size %d!\n", (int)size);
        }*/
    };
}

#endif
