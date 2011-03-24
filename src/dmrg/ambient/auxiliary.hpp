#include "ambient/auxiliary.h"
namespace ambient{

    template<typename T>
    operation_stack<T>::operation_stack() : write_iterator(0), read_iterator(0), alt_read_iterator(0), length(0) {
        this->content = (T*)malloc(sizeof(T)*STACK_CONTENT_RESERVATION);
        this->reserved = STACK_CONTENT_RESERVATION;
    }

    template<typename T>
    void operation_stack<T>::push_back(T element){
        this->content[this->write_iterator++] = element;
        this->length++;
        if(this->write_iterator == this->reserved){
            this->reserved += STACK_CONTENT_RESERVATION;
            this->content = (T*)realloc(this->content, sizeof(T)*this->reserved);
        }
    }

    template<typename T>
    bool operation_stack<T>::end_reached(){
        if(this->read_iterator == this->length){
            this->read_iterator = this->write_iterator = 0;
            return true;
        }
        return false;
    }

    template<typename T>
    bool operation_stack<T>::alt_end_reached(){
        if(this->alt_read_iterator == this->length){
            this->alt_read_iterator = 0;
            return true;
        }
        return false;
    }

    template<typename T>
    void operation_stack<T>::clean(){
        this->length = 0;
    }

    template<typename T>
    T* operation_stack<T>::pick(){
        return &this->content[this->read_iterator++];
    }

    template<typename T>
    T* operation_stack<T>::alt_pick(){
        return &this->content[this->alt_read_iterator++];
    }
}
