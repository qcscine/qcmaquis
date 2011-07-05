#ifndef __VLI_BANK_H
#define __VLI_BANK_H

#include "align.h"

namespace bank {

#define MAXMEM_CPU 1024    
#define MAXMEM_GPU 1024    // to do change by  engine
    
template<class VLI>
class membank{
public:    
    typedef typename VLI::value_type value_type; 
    typedef typename VLI::size_type size_type;
    membank(){};
    virtual ~membank(){};
    virtual void* malloc(size_type buffer)=0;     
    virtual void  free(value_type* mem_buffer)=0;
};

template<class VLI>
class membankcpu: public membank<VLI>{
public:
    typedef typename VLI::value_type value_type; 
    typedef typename VLI::size_type size_type;
    
    membankcpu():offset_(0),previouspos_(-1),ressource_bank_(MAXMEM_CPU*sizeof(value_type)){
        mem_ = (value_type*)align::aligned_malloc(ressource_bank_); // cast void to the type of the bas
    }
    
    ~membankcpu(){
        align::aligned_free(this->mem_);   
    }
    /** take a buffer into the bank e.g. a polynomial**/
    void* malloc(size_type buffer){
        value_type* mem_offset = (mem_+offset_);
        offset_ += (buffer/static_cast<size_type>(sizeof(value_type)));
        ressource_bank_ -= buffer;
        if (ressource_bank_ < 0){
            //except EX("CPU Bankrupt, GAME OVER ! ^^' "); TO DO !
        }
        if(previouspos_ != -1){ 
            offset_ = previouspos_;
            previouspos_ = -1;
        }
        return ((void*)(mem_offset));
    }
    
    void free(value_type* mem_buffer){
        previouspos_ = offset_;
        offset_ = ((size_type)mem_buffer-(size_type)mem_)/static_cast<size_type>(sizeof(value_type));
    }
    
    
private:
    value_type* mem_;
    size_type offset_, previouspos_;
    long int ressource_bank_;
};
    
} // end namespace

#endif