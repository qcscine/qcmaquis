#define COLLECTOR_STR_RESERVE 65536
#define COLLECTOR_RAW_RESERVE 1024

namespace ambient{

    using ambient::models::velvet::history;
    using ambient::models::velvet::revision;

    inline collector::collector(){
        this->str.reserve(COLLECTOR_STR_RESERVE);
        this->raw.reserve(COLLECTOR_RAW_RESERVE);
    }

    inline void collector::push_back(void* o){
        this->raw.push_back(o);
    }

    inline void collector::push_back(history* o){
        this->str.push_back(o);
    }

    inline void collector::delete_ptr::operator()( history* element ) const {
        int size = element->content.size();
        for(int i = 0; i < size; i++) ambient::controller.free(*element->content[i]);
        for(int i = 0; i < size; i++) ambient::static_memory::free<revision>(element->content[i]); 
        delete element;
    }

    inline void collector::delete_ptr::operator()( void* element ) const {
        ambient::static_memory::free<FUTURE_SIZE>(element);
    } 

    inline void collector::clear(){
        std::for_each( str.begin(), str.end(), delete_ptr());
        std::for_each( raw.begin(), raw.end(), delete_ptr());
        str.clear();
        raw.clear();
    }

}

