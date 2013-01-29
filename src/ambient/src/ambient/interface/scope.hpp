#ifndef AMBIENT_INTERFACE_SCOPE
#define AMBIENT_INTERFACE_SCOPE

namespace ambient { namespace scope {

    class single { 
    public:
        //single(){
        //    printf("Created Ambient individual scope!\n");
        //}
        //single(size_t rank){
        //    printf("Created Ambient individual scope from rank %d!\n", (int)rank);
        //}
    };

    class shared { 
    public:
        shared(){
            printf("Created Ambient shared scope!\n");
        }
        shared(size_t size){
            printf("Created Ambient shared scope of size %d!\n", (int)size);
        }
    };

    class common { 
    public:
        common(){
            printf("Falling back to the common/world scope!\n");
        }
    };

} }

#endif
