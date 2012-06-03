#ifndef AMBIENT_UTILS_SINGLETON
#define AMBIENT_UTILS_SINGLETON

namespace ambient {

    template <typename T>
    class singleton {
    public:
        static T& instance();
        inline virtual ~singleton(){ }          // execute derived destructors
    protected:
        inline singleton(){}                    // only for derived classes
    private:
        singleton(singleton const&);            // copy constructor is private
        singleton& operator=(singleton const&); // assignment operator is private
    };

    template <typename T>
    inline T& singleton<T>::instance(){
        static T singleton;                     // not thread-safe
        return singleton;
    }
}

#endif
