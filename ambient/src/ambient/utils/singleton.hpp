#ifndef AMBIENT_SINGLETON_HPP
#define AMBIENT_SINGLETON_HPP

namespace ambient {

    template <typename T>
    class singleton {
    public:
        static T& instance();
        virtual ~singleton(){ }                 // execute derived destructors
    protected:
        singleton(){}                           // only for derived classes
    private:
        singleton(singleton const&);            // copy constructor is private
        singleton& operator=(singleton const&); // assignment operator is private
    };

    template <typename T>
    T& singleton<T>::instance(){
        static T singleton;                     // not thread-safe
        return singleton;
    }

}

#endif
