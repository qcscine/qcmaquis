#include "ambient/utils/hashmap.h"

namespace ambient{

    hashmap::hashmap()
    : content(HASHMAP_PARTITION_SIZE)
    {
    }
    
    size_t hashmap::insert(size_t* hash, size_t hash_len, void* value, size_t shift){
        size_t hash_w = hash[0] >> shift;
        if(hash_w >= HASHMAP_PARTITION_SIZE){
            size_t hash_cut = (size_t)(unsigned char)hash_w; // first log2 of HASHMAP_PARTITION_SIZE bits
            if(this->content[hash_cut].first == NULL){
                this->content[hash_cut].first = new hashmap();
            }
            return this->content[hash_cut].first->insert(hash, hash_len, value, shift+HASHMAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            if(this->content[hash_w].first == NULL){
                this->content[hash_w].first = new hashmap();
            }
            return this->content[hash_w].first->insert(&hash[1], hash_len-1, value); 
        }
    
        if(this->content[hash_w].second == NULL){
            this->content[hash_w].second = new std::vector<void*>();
        }
        std::vector<void*>* data = this->content[hash_w].second;
        if(data->capacity() == data->size())
            data->reserve(data->size()+HASHMAP_VECTOR_RESERVATION);
        data->resize(data->size()+1);
        
        (*data)[data->size()-1] = value;
        return data->size(); // 1-based id
    }
    
    void* hashmap::find(size_t* hash, size_t hash_len, size_t id, size_t shift) const {
        size_t hash_w = hash[0] >> shift;
        if(hash_w >= HASHMAP_PARTITION_SIZE){
            size_t hash_cut = (size_t)(unsigned char)hash_w;
            return this->content[hash_cut].first->find(hash, hash_len, id, shift+HASHMAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            return this->content[hash_w].first->find(&hash[1], hash_len-1, id); 
        }
    
        return (*this->content[hash_w].second)[id-1];
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
