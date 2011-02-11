#include "ambient/auxiliary.h"

namespace ambient{

    template <class T>
    hash_map<T>::hash_map():content(HASH_MAP_PARTITION_SIZE){ }
    
    template <class T>
    void hash_map<T>::insert(unsigned int* hash, unsigned int hash_len, unsigned int id, T value, int shift)
    {
        unsigned int hash_w = hash[0] >> shift;
        if(hash_w >= HASH_MAP_PARTITION_SIZE){
            unsigned int hash_cut = (unsigned int)(unsigned char)hash_w; // first log2 of HASH_MAP_PARTITION_SIZE bits
            if(this->content[hash_cut].first == NULL){
                this->content[hash_cut].first = new hash_map();
            }
            return this->content[hash_cut].first->insert(hash, hash_len, id, value, shift+HASH_MAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            if(this->content[hash_w].first == NULL){
                this->content[hash_w].first = new hash_map();
            }
            return this->content[hash_w].first->insert(&hash[1], hash_len-1, id, value); 
        }
    
        if(this->content[hash_w].second == NULL){
            this->content[hash_w].second = new std::vector<T>(HASH_MAP_VECTOR_RESERVATION);
        }
        std::vector<T>* data = this->content[hash_w].second;
        if(data->capacity() < id){
            data->reserve(id+HASH_MAP_VECTOR_RESERVATION);
            data->resize(id);
        }else if(data->size() < id){
            data->resize(id);
        }
        (*data)[id] = value;
    }
    
    template <class T>
    T hash_map<T>::find(unsigned int* hash, unsigned int hash_len, unsigned int id, int shift)
    {
        unsigned int hash_w = hash[0] >> shift;
        if(hash_w >= HASH_MAP_PARTITION_SIZE){
            unsigned int hash_cut = (unsigned int)(unsigned char)hash_w;
            return this->content[hash_cut].first->find(hash, hash_len, id, shift+HASH_MAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            return this->content[hash_w].first->find(&hash[1], hash_len-1, id); 
        }
    
        return (*this->content[hash_w].second)[id];
    }

    hash_map<core::coherency_table> void_pt_map;

}
