#include "ambient/auxiliary.h"
#include "ambient/core/operation/operation.h"
#include <stdlib.h>

namespace ambient{

    access_marker::access_marker& access_marker::instance()
    {
        static access_marker* singleton = NULL;
        if(!singleton) singleton = new access_marker();
        return *singleton;
    }
    access_marker::access_marker():write_only(false),nesting(0){ }
    
    void access_marker::write_only_start_mark(){
        if(this->nesting++ == 0) this->write_only = true;
    }
    void access_marker::write_only_stop_mark(){
        if(--this->nesting == 0) this->write_only = false;
    }
    bool access_marker::write_only_marked(){ return this->write_only; }
    
    hash_map::hash_map& hash_map::instance()
    {
        static hash_map* singleton = NULL;
        if(!singleton) singleton = new hash_map();
        return *singleton;
    }
    hash_map::hash_map():content(HASH_MAP_PARTITION_SIZE){ }
    
    unsigned int hash_map::insert(unsigned int* hash, unsigned int hash_len, core::layout_table* value, int shift)
    {
        unsigned int hash_w = hash[0] >> shift;
        if(hash_w >= HASH_MAP_PARTITION_SIZE){
            unsigned int hash_cut = (unsigned int)(unsigned char)hash_w; // first log2 of HASH_MAP_PARTITION_SIZE bits
            if(this->content[hash_cut].first == NULL){
                this->content[hash_cut].first = new hash_map();
            }
            return this->content[hash_cut].first->insert(hash, hash_len, value, shift+HASH_MAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            if(this->content[hash_w].first == NULL){
                this->content[hash_w].first = new hash_map();
            }
            return this->content[hash_w].first->insert(&hash[1], hash_len-1, value); 
        }
    
        if(this->content[hash_w].second == NULL){
            this->content[hash_w].second = new std::vector<core::layout_table*>();
        }
        std::vector<core::layout_table*>* data = this->content[hash_w].second;
        if(data->capacity() == data->size())
            data->reserve(data->size()+HASH_MAP_VECTOR_RESERVATION);
        data->resize(data->size()+1);
        
        (*data)[data->size()-1] = value;
        return data->size(); // 1-based id
    }
    
    core::layout_table* hash_map::find(unsigned int* hash, unsigned int hash_len, unsigned int id, int shift) const
    {
        unsigned int hash_w = hash[0] >> shift;
        if(hash_w >= HASH_MAP_PARTITION_SIZE){
            unsigned int hash_cut = (unsigned int)(unsigned char)hash_w;
            return this->content[hash_cut].first->find(hash, hash_len, id, shift+HASH_MAP_PARTITION_BIT_SIZE); 
        }else if(hash_len > 1){
            return this->content[hash_w].first->find(&hash[1], hash_len-1, id); 
        }
    
        return (*this->content[hash_w].second)[id-1];
    }

    delegate::delegate() : length(0) {
        this->handlers = (core::operation**)malloc(sizeof(core::operation*)*STACK_CONTENT_RESERVATION);
        this->reserved = STACK_CONTENT_RESERVATION;
    }

    void delegate::operator+=(core::operation* handler){
        this->handlers[this->length++] = handler;
        if(this->length == this->reserved){
            this->reserved += STACK_CONTENT_RESERVATION;
            this->handlers = (core::operation**)realloc(this->handlers, sizeof(core::operation*)*this->reserved);
        }
    }

    void delegate::operator()(){
        for(int i=0; i < this->length; i++)
            this->handlers[i]->invoke();
    }
}
