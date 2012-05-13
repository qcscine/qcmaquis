#include "ambient/utils/hashmap.h"
#define HASHMAP_VECTOR_RESERVATION 1.5

namespace ambient{

    hashmap::hashmap()
    : content(HASHMAP_PARTITION_SIZE)
    {
    }
    
    size_t hashmap::insert(void* value){
        size_t size = this->content.size();
        if(this->content.capacity() == size)
            this->content.reserve(size*HASHMAP_VECTOR_RESERVATION);
        this->content.resize(size+1);
        
        this->content[size] = value;
        return size+1; // 1-based id
    }
    
    fhashmap::fhashmap()
    : content(HASHMAP_PARTITION_SIZE)
    { 
    }
    
    void** fhashmap::get(size_t* hash, size_t hash_len, size_t shift) const {
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

}
