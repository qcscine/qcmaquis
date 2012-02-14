#ifndef AMBIENT_INTERFACE_CONTROLLER_H
#define AMBIENT_INTERFACE_CONTROLLER_H
#include "ambient/models/imodel.h"
#include "ambient/channels/ichannel.h"

namespace ambient { namespace controllers { 

    class icontroller {
    public:
        virtual void acquire(channels::ichannel* channel) = 0;
        virtual models::imodel::layout::entry& fetch_block(models::imodel::revision& r, size_t i, size_t j) = 0;
        virtual models::imodel::layout::entry& ifetch_block(models::imodel::revision& r, size_t i, size_t j) = 0;
        virtual void push(models::imodel::modifier* op) = 0;
        virtual void atomic_complete() = 0;
        virtual void flush() = 0;
    };
    
    void forward_block(channels::ichannel::packet&);
    void accept_block(channels::ichannel::packet&);
} }

namespace ambient {
    extern controllers::icontroller& controller;
}

#endif
