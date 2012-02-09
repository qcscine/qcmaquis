#include "ambient/utils/touchstack.h"

namespace ambient{

    template<typename T>
    touchstack<T>::touchstack() 
    : write_iterator(0), read_iterator(0), alt_read_iterator(0), length(0) 
    {
        this->content = (T*)malloc(sizeof(T)*STACK_CONTENT_RESERVATION);
        this->reserved = STACK_CONTENT_RESERVATION;
    }

    template<typename T>
    touchstack<T>::~touchstack(){
        free(this->content);
    }

    template<typename T>
    void touchstack<T>::push_back(T element){
        this->content[this->write_iterator++] = element;
        this->length++;
        if(this->write_iterator == this->reserved){
            this->reserved += STACK_CONTENT_RESERVATION;
            this->content = (T*)realloc(this->content, sizeof(T)*this->reserved);
        }
    }

    template<typename T>
    bool touchstack<T>::end_reached(){
        if(this->read_iterator == this->length){
            this->read_iterator = this->write_iterator = 0;
            return true;
        }
        return false;
    }

    template<typename T>
    void touchstack<T>::reset(){
        this->read_iterator = this->write_iterator = 0;
    }

    template<typename T>
    void touchstack<T>::alt_reset(){
        this->alt_read_iterator = 0;
    }

    template<typename T>
    bool touchstack<T>::alt_end_reached(){
        if(this->alt_read_iterator == this->length){
            this->alt_read_iterator = 0;
            return true;
        }
        return false;
    }

    template<typename T>
    void touchstack<T>::clean(){
        this->length = 0;
        this->read_iterator = 0;
        this->alt_read_iterator = 0;
        this->write_iterator = 0;
    }

    template<typename T>
    bool touchstack<T>::empty(){
        return (this->length == 0);
    }

    template<typename T>
    void touchstack<T>::sort(){ // insertion sort
        for(int i = 1; i < this->length; i++){
            T value = this->content[i];
            int j = i - 1;
            bool done = false;
            do {
                if(this->content[j]->get_weight() < value->get_weight()){
                    this->content[j + 1] = this->content[j];
                    if(--j < 0) done = true;
                }else
                    done = true;
            } while(!done);
            this->content[j + 1] = value;
        }
    }

    template<typename T>
    T touchstack<T>::pick(){
        return this->content[this->read_iterator++];
    }

    template<typename T>
    T touchstack<T>::alt_pick(){
        return this->content[this->alt_read_iterator++];
    }

    template<typename T>
    T touchstack<T>::back(){
        return this->content[this->length-(++this->read_iterator)];
    }

} // namespace ambient
