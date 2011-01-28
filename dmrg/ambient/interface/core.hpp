    template <typename T> p_profile::p_profile(const T* ptr){ 
        p_profile_model(this, ptr); 
        for(int j = 0; j < this->dim.x / scheduler::instance().group_dim().x; j++) // fortran matrix style ,)
            for(int i = 0; i < this->dim.y / scheduler::instance().group_dim().y; i++)
                for(int k = 0; k < this->dim.z / scheduler::instance().group_dim().z; k++)
                    skeleton.push_back(new workgroup(this, i, j, k));
    };

    template <typename RT, typename FC, typename FL, class T1, class T2> 
    const RT push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2){
	zout << "Calling the template that takes 2 args\n";
	p_profile* handle = new p_profile((const RT*)NULL);
	RT out(handle);
	push(c_kernel, l_kernel, arg1, arg2, out);
	return RT(handle);
    }

    template <typename FC, typename FL, class T1, class T2, class T3>
    void push(FC c_kernel, FL l_kernel, const T1& arg1, const T2& arg2, const T3& arg3){
	l_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3)); 
	c_kernel(get_profile(arg1), get_profile(arg2), get_profile(arg3));
    }

    template <typename L, typename R>
    void pin(const L& lhs, const R& rhs){
//        new p_action('=', &lhs, &rhs);
    }

    template <typename L, typename R>
    void copy(const L& lhs, const R& rhs){
//        new p_action('|', &lhs, &rhs);
    }
