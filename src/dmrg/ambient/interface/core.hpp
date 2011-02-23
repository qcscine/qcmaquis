namespace ambient {
using namespace blas;

    #include "ambient/interface/push.pp.hpp"

    template <typename T> void_pt::void_pt(const T* ptr) 
    : p_profile()
    {
        void_pt_model(this, ptr);
        this->regroup(); 
    };

    template <typename L, typename R>
    void pin(L& proxy_object, const R& real_object){
        void_pt* proxy = get_profile(proxy_object);
        void_pt* real  = get_profile(real_object);
        proxy->profile   = real->profile;
        proxy->set_dim(real->get_dim());
        real->imitate(proxy); // copy proxy settings to the actual profile
    }
}
