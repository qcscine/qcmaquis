#ifndef AMBIENT_HASHMAP_H
#define AMBIENT_HASHMAP_H

#define HASHMAP_PARTITION_BIT_SIZE 8
#define HASHMAP_PARTITION_SIZE 256 // 2^HAH_MAP_PARTITION_BIT_SIZE
#define HASHMAP_VECTOR_RESERVATION 10

#include <stdlib.h>
#include <vector>

namespace ambient{

    class hashmap {
    public:
        hashmap();
        size_t insert(size_t* hash, size_t hash_len, void* value, size_t shift = 0);
        void* find(size_t* hash, size_t hash_len, size_t id, size_t shift = 0) const;
    private:
        std::vector< std::pair<hashmap*,std::vector<void*>* > > content;
    };

    class fhashmap {
    public:
        fhashmap();
        void** get(size_t* hash, size_t hash_len, size_t shift = 0) const;
        mutable std::vector< std::pair<fhashmap*, void*> > content;
    };

} // namespace ambient
#endif
