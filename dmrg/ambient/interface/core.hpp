namespace ambient { using namespace blas;

    template <typename FL, typename FC, class T0, class T1, class T2>
    void push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1, T2& arg2){
        if(get_profile(arg0)->proxy) pin(arg0, arg2);
        if(get_profile(arg1)->proxy) pin(arg1, arg2);
        ambient::engine.push(new core::operation(l_kernel, &arg0, &arg1, &arg2),
                             new core::operation(c_kernel, &arg0, &arg1, &arg2));
    }
    template <typename ST, typename FL, typename FC, class T0, class T1>
    ST push(FL l_kernel, FC c_kernel, T0& arg0, T1& arg1){
        void_pt* handle = new void_pt((ST*)NULL);
        ST out(handle);
        push(l_kernel, c_kernel, arg0, arg1, out);
        return out;
    }
    #include "ambient/interface/push.pp.hpp" // all other variants

    template <typename L, typename R>
    void pin(L& proxy_object, const R& real_object){
        void_pt* proxy = get_profile(proxy_object);
        void_pt* real  = get_profile(real_object);
        proxy->profile   = real->profile;
        proxy->set_dim(real->get_dim());
        real->imitate(proxy); // copy proxy settings to the actual profile
    }
}
