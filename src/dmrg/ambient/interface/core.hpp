// nested inside ambient.hpp in ambient namespace
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

template <typename T> 
void_pt& breakdown(T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.profile);
    *profile_ptr = (void_pt*)(*profile_ptr)->dereference();
    (*profile_ptr)->inconstant();
    return **profile_ptr;
}

template <typename T> 
void_pt& breakdown(const T& obj){
    void_pt** profile_ptr = const_cast<void_pt**>(&obj.profile);
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
void plus_reduce(workgroup* grp, void* update){
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
