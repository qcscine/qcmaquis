    template <typename T> p_profile::p_profile(const T* ptr) : specific(false), profile(this) { 
        p_profile_model(this, ptr);
        this->regroup(); 
    };

    template <typename ST, typename FC, typename FL, class T1, class T2> 
    const ST push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2){
	zout << "Calling the template that takes 2 args\n";
	p_profile* handle = new p_profile((const ST*)NULL);
	ST out(handle);
	push(c_kernel, l_kernel, arg1, arg2, out);
	return out;
    }

    template <typename FC, typename FL, class T1, class T2, class ST>
    void push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2, const ST& structuring_arg){
        if(get_profile(arg1)->proxy) pin(const_cast<T1&>(arg1), structuring_arg);
        if(get_profile(arg2)->proxy) pin(const_cast<T2&>(arg2), structuring_arg);

        core::operation* l_operation = new core::operation(l_kernel, get_profile(arg1), get_profile(arg2), get_profile(structuring_arg)); 
        l_operation->perform();
        core::operation* c_operation = new core::operation(c_kernel, get_profile(arg1), get_profile(arg2), get_profile(structuring_arg)); 
        c_operation->perform();
    }

    template <typename L, typename R>
    void pin(L& proxy_object, const R& real_object){
        zout << "pinning arguments \n";
        p_profile* proxy = get_profile(proxy_object);
        p_profile* real  = get_profile(real_object);
        proxy->profile   = real->profile;
        proxy->dim       = real->dim;
// copy proxy settings to the actual profile
        real->specific   = proxy->specific;
        real->dim_distr  = proxy->dim_distr;
        real->dim_group  = proxy->dim_group;
        real->dim_item   = proxy->dim_item;
        real->dim_gpu    = proxy->dim_gpu;
        real->regroup();
    }

    template <typename L, typename R>
    void copy(const L& lhs, const R& rhs){
//        new p_action('|', &lhs, &rhs);
    }
