#ifndef AMBIENT_UTILS_HASHMAP
#define AMBIENT_UTILS_HASHMAP
#define HASHMAP_PARTITION_BIT_SIZE 8
#define HASHMAP_PARTITION_SIZE 256 // 2^HAH_MAP_PARTITION_BIT_SIZE
#define HASHMAP_VECTOR_RESERVATION 10
#include <stdlib.h>
#include <vector>

namespace ambient{

    class hashmap {
    public:
        hashmap();
        size_t insert(void* value);
        inline void* find(size_t id) const { return this->content[id-1]; };
    private:
        std::vector<void*> content;
    };

    class fhashmap {
    public:
        fhashmap();
        void** get(size_t* hash, size_t hash_len, size_t shift = 0) const;
        mutable std::vector< std::pair<fhashmap*, void*> > content;
    };

}

#endif
