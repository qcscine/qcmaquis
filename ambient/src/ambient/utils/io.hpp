#ifndef AMBIENT_IO
#define AMBIENT_IO

#include <boost/tuple/tuple.hpp>

namespace ambient {

    template<typename T> class future;
    bool verbose();

}

namespace ambient { namespace utils {

    class fstream {
    public: 
         fstream(){ 
             vec_of_list_tiles_ = std::vector<std::list< boost::tuple<void*, size_t, size_t> > >(__cilkrts_get_total_workers()); // idea : each cilk thread saves in a seperate list, I avoid mutex
         };

         void sort(std::list<boost::tuple<void*, size_t, size_t> >& full_list_tiles){
             typename std::list< boost::tuple<void*, size_t, size_t> >::iterator it_list = full_list_tiles.begin();
             typename std::vector< std::list< boost::tuple<void*, size_t, size_t> > >::iterator it_vec = this->vec_of_list_tiles_.begin();
             for(; it_vec != this->vec_of_list_tiles_.end(); ++it_vec)
                 full_list_tiles.splice(it_list, *it_vec); // with splice I clean the vec of list 

             full_list_tiles.sort([](const boost::tuple<void*, size_t, size_t>& lhs, boost::tuple<void*, size_t, size_t>& rhs)
                                    { return boost::get<1>(lhs) < boost::get<1>(rhs); } ); // pair are sorted in function of the 2nd argument of the tuple, c++11 lambda features
         }

         void save(void* buffer, size_t tag, size_t size_buffer){ // call by ambient kernel
             build_list_write(buffer, tag, size_buffer);
         }

         void load(void* buffer, size_t tag, size_t size_buffer){ // call by ambient kernel
             build_list_read(buffer, tag, size_buffer);
         }

         void build_list_write(void* buffer, size_t tag, size_t size_buffer){
             void* value = malloc(size_buffer);
             memcpy(value, buffer, size_buffer); // I can not create a temporary matrix by an independant thread without crash so cpy ..., to change one day
             this->vec_of_list_tiles_[__cilkrts_get_worker_number()].push_back(boost::make_tuple(value, tag, size_buffer));
         }

         void build_list_read(void* buffer, size_t tag, size_t size_buffer){
             this->vec_of_list_tiles_[__cilkrts_get_worker_number()].push_back(boost::make_tuple(buffer, tag, size_buffer));
         }

         void get_list(std::list< boost::tuple<void*, size_t, size_t> >& full_list_tiles){ 
             sort(full_list_tiles);
         }   

    private:
         std::vector<std::list< boost::tuple<void*,size_t, size_t> > > vec_of_list_tiles_; // we are working with memory, so void*, if desagree -> Alex
    };

    class ostream {
    public:
        std::fstream nullio;
        ostream() : nullio("/dev/null") { }

        template<class T>
        ostream& operator<<(future<T> const & obj){
            std::cout << obj.calc_value();
            return *this;
        }

        template<class T>
        ostream& operator<<(T const & obj){
            std::cout << obj;
            return *this;
        }

        ostream& operator<<(std::ostream& (*pf)(std::ostream&)){
            std::cout << pf;
            return *this;
        }

        void precision(int p){
            std::cout.precision(p);
        }

        void flush(){
            std::cout.flush();
        }
    };

    class mpostream {
    public:
        std::fstream nullio;
        mpostream() : nullio("/dev/null") { }

        template<class T>
        mpostream& operator<<(future<T> const & obj){
            if(ambient::verbose()) std::cout << obj.calc_value();
            else nullio << obj.calc_value();
            return *this;
        }

        template<class T>
        mpostream& operator<<(T const & obj){
            if(ambient::verbose()) std::cout << obj;
            return *this;
        }

        mpostream& operator<<(std::ostream& (*pf)(std::ostream&)){
            if(ambient::verbose()) std::cout << pf;
            return *this;
        }

        void precision(int p){
            if(ambient::verbose()) std::cout.precision(p);
        }

        void flush(){
            if(ambient::verbose()) std::cout.flush();
        }
    };

} }

namespace ambient {
    extern utils::fstream fout;
    extern utils::mpostream cout;
    extern utils::mpostream cerr;
}

#endif
