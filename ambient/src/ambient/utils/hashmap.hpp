#ifndef AMBIENT_UTILS_HASHMAP
#define AMBIENT_UTILS_HASHMAP
#define HASHMAP_PARTITION_BIT_SIZE 8
#define HASHMAP_PARTITION_SIZE 65536 // 2^HAH_MAP_PARTITION_BIT_SIZE
#define HASHMAP_VECTOR_RESERVATION 2

namespace ambient{

    class hashmap {
    public:
        hashmap():content(HASHMAP_PARTITION_SIZE){}
        inline size_t insert(void* value){
            size_t size = this->content.size();
            if(this->content.capacity() == size)
                this->content.reserve(size*HASHMAP_VECTOR_RESERVATION);
            this->content.push_back(value);
            return size;
        }
        inline void* find(size_t id) const { return this->content[id]; };
    private:
        std::vector<void*> content;
    };

    class fhashmap {
    public:
        fhashmap():content(HASHMAP_PARTITION_SIZE){}
        inline void** get(size_t* hash, size_t hash_len, size_t shift = 0) const {
            size_t hash_w = hash[0] >> shift;
            if(hash_w >= HASHMAP_PARTITION_SIZE){
                size_t hash_cut = (size_t)(unsigned char)hash_w; // first log2 of HASHMAP_PARTITION_SIZE bits
                if(this->content[hash_cut].first == NULL){
                    this->content[hash_cut].first = new fhashmap();
                }
                return this->content[hash_cut].first->get(hash, hash_len, shift+HASHMAP_PARTITION_BIT_SIZE); 
            }else if(hash_len > 1){
                if(this->content[hash_w].first == NULL){
                    this->content[hash_w].first = new fhashmap();
                }
                return this->content[hash_w].first->get(&hash[1], hash_len-1); 
            }
            return &this->content[hash_w].second;
        }
        mutable std::vector< std::pair<fhashmap*, void*> > content;
    };

}

#endif
