namespace ambient {
using namespace blas;

    template <typename T> void_pt::void_pt(const T* ptr) 
    : p_profile()
    {
        void_pt_model(this, ptr);
        this->regroup(); 
    };

    template <typename FL, typename FC, class T1, class T2, class T3>
    void push(FL l_kernel, FC c_kernel, T1& arg1, T2& arg2, T3& arg3){
        if(get_profile(arg1)->proxy) pin(arg1, arg3);
        if(get_profile(arg2)->proxy) pin(arg2, arg3);

        ambient::engine.push(new core::operation(l_kernel, &arg1, &arg2, &arg3), 
                             new core::operation(c_kernel, &arg1, &arg2, &arg3)); 
    }

    template <typename ST, typename FL, typename FC, class T1, class T2> 
    ST push(FL l_kernel, FC c_kernel, T1& arg1, T2& arg2){
	void_pt* handle = new void_pt((ST*)NULL);
	ST out(handle);
	push(l_kernel, c_kernel, arg1, arg2, out);
	return out;
    }

    template <typename L, typename R>
    void pin(L& proxy_object, const R& real_object){
        void_pt* proxy = get_profile(proxy_object);
        void_pt* real  = get_profile(real_object);
        proxy->profile   = real->profile;
        proxy->set_dim(real->get_dim());
        real->imitate(proxy); // copy proxy settings to the actual profile
    }
}
