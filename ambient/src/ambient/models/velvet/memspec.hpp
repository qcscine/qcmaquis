#define BOUND 30

namespace ambient { namespace models { namespace velvet {

    inline memspec::memspec(size_t size, dim2 block, dim2 dim)
    : clean(false), size(size), block(block), grid(__a_ceil(dim.x/block.x), __a_ceil(dim.y/block.y)), dim(dim)
    {
    }

    inline void* memspec::alloc() const {    
        return malloc(size + BOUND);
    }

    inline void* memspec::calloc() const {    
        void* memory = malloc(size + BOUND);
        memset((char*)memory + BOUND, 0, size);
        return memory;
    }

    inline size_t memspec::get_bound() const {
        return BOUND; // ambient::channel.get_block_packet_type(size).get_bound(9);
    }
 
    void* memspec::get_packet_t() const {
        return &ambient::channel.get_block_packet_type(size);
    }

} } }

#undef BOUND
