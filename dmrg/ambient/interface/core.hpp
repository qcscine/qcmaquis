// nested inside ambient.hpp in ambient namespace
using namespace blas;

#include "ambient/core/operation/operation.pp.hpp"
#define MAX_NUM_CHAR_LEN 10
#define scope_select(...) scope_select(std::string(std::string() + __VA_ARGS__).c_str());

class void_pt: public p_profile 
{ 
public: 
    template <typename T> void_pt(const T* ptr) : p_profile()
    { breakdown_model(this, ptr); }
private: 
    ~void_pt(){ };  // making user unable to delete the profile
    friend class T; // the container can delete its profile
    template<class T> friend inline void boost::checked_delete(T * x); // so as boost >_<
};

template<typename T>
void null_deleter(T* object){ }

template<typename T>
void poke_deleter(T* object){
//    printf("Poking... %d - %d ---- proxy: %d\n", object->use_count, ambient::rank(), object->breakdown()->id);
    if(--object->use_count) return;
    delete object;
}

template<typename T, policy P>
class livelong
{
public:
    livelong():p(P),use_count(0){
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica()) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    livelong(void_pt* p): p(P),use_count(0),profile(p)
    {
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(p)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    template<typename A1>
    livelong(A1 a1):p(P),use_count(0){
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    template<typename A1, typename A2>
    livelong(A1 a1, A2 a2):p(P),use_count(0){
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1, a2)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    template<typename A1, typename A2, typename A3>
    livelong(A1 a1, A2 a2, A3 a3):p(P),use_count(0){
        if(P == ANY) handle.reset( self = (T*)(new typename T::replica(a1, a2, a3)) );
        else if(P == MANUAL) handle.reset( self = (T*)this, null_deleter<T> );
        else if(P == REPLICA || P == WEAK) self = (T*)this;
        thyself = self;
    }
    boost::shared_ptr<T> get_handle() const {
        if(this->p == ANY) return this->handle;
        else if(this->p == MANUAL) return this->handle;
        else if(this->p == WEAK){ ((livelong*)this)->use_count++; return boost::shared_ptr<T>(this->self, &poke_deleter<T>); } // for one-touch only
        else if(this->p == REPLICA) return boost::shared_ptr<T>(this->self);
        return boost::shared_ptr<T>();
    }
    void_pt*& breakdown() const {
        return self->profile;
    }
    void_pt*& set_breakdown() {
        return (self->profile = new void_pt(self));
    }
public:
    T* self;
    T* thyself;
    size_t use_count;
private:
    policy p; // the same as P (avoiding casting collisions)
    void_pt* profile;
    boost::shared_ptr<T> handle;
};

template<typename T>
boost::shared_ptr<T> get_handle(T& a){
    return boost::shared_ptr<T>(new T(a));
}

template<typename T>
boost::shared_ptr<T> get_handle(const T& a){
    return boost::shared_ptr<T>(new T(a));
}

template<typename T, policy P>
boost::shared_ptr< p_dense_matrix<T> > get_handle(p_dense_matrix<T,P>& a){
    return boost::static_pointer_cast< p_dense_matrix<T> >(boost::static_pointer_cast<void>(a.get_handle()));
}

template<typename T, policy P>
boost::shared_ptr< p_dense_matrix<T> > get_handle(const p_dense_matrix<T,P>& a){
    return boost::static_pointer_cast< p_dense_matrix<T> >(boost::static_pointer_cast<void>(a.get_handle()));
}

template <typename T> 
void_pt& breakdown(T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.breakdown());
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->inconstant();
    return **profile_ptr;
}

template <typename T> 
void_pt& breakdown(const T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.breakdown());
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->constant();
    return **profile_ptr;
}

// aliases for computational kernels
template <typename T>
void_pt& current(T& obj){
    return breakdown(obj);
}

template <typename T>
void plus_reduce(memblock* grp, void* update){
    assert(false); // only partially specialized
}

template <char R, typename T>
void_pt& reduced(T& obj){
    if(breakdown(obj).associated_proxy == NULL){
        if(R == '+') breakdown(obj).associate_proxy(new void_pt((T*)NULL), plus_reduce<T>);
        breakdown_proxy_model((void_pt*)breakdown(obj).associated_proxy, &(breakdown(obj)), &obj);
    }
    return *((void_pt*)breakdown(obj).associated_proxy);
}

std::string & operator+(std::string & lhs, double rhs){
    char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
    if(rhs - (int)rhs) sprintf(rhs_str,"%.1f",rhs);
    else sprintf(rhs_str,"%d",(int)rhs);
    lhs += rhs_str;
    free(rhs_str);
    return lhs;
}
std::string & operator+(std::string & lhs, int rhs){
    char* rhs_str = (char*)malloc(sizeof(char)*MAX_NUM_CHAR_LEN);
    sprintf(rhs_str,"%d", rhs);
    lhs += rhs_str;
    free(rhs_str);
    return lhs;
}
std::string & operator+(const std::string & lhs, double rhs){
    return const_cast<std::string&>(lhs)+rhs;
}
std::string & operator+(const std::string & lhs, int rhs){
    return const_cast<std::string&>(lhs)+rhs;
}
std::string & operator+(std::string & lhs, const char* rhs){
    return lhs += rhs;
}
std::string & operator+(const std::string & lhs, const char* rhs){
    return const_cast<std::string&>(lhs) += rhs;
}
std::string & operator+(std::string & lhs, std::pair<unsigned int*,size_t> rhs){
    lhs + (int)(*rhs.first);
    lhs += "-";
    lhs + (int)(rhs.second);
    return lhs;
}
std::string & operator+(const std::string & lhs, std::pair<unsigned int*,size_t> rhs){
    const_cast<std::string&>(lhs) + (int)(*rhs.first);
    const_cast<std::string&>(lhs) += "-";
    return const_cast<std::string&>(lhs) + (int)rhs.second;
}
