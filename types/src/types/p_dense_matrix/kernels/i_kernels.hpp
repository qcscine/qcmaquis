#ifndef __MAQUIS_TYPES_I_KERNELS_HPP__
#define __MAQUIS_TYPES_I_KERNELS_HPP__
namespace ambient {

    template<typename T> void randomize(T* ad){ *ad = drand48(); }
    template<typename T> void randomize(std::complex<T>* ad){
        (*ad).real() = drand48();
        (*ad).imag() = drand48();
    }

    template<typename T>
    void random_i(models::v_model::object& a){
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t m = get_mem_dim(a).y;
        size_t n = get_mem_dim(a).x;
        size_t ld = m;
        size_t sd = n;
        T* ad = current(a)(i,j);
    
        if(i == (get_grid_dim(a).y-1) && get_dim(a).y%get_mem_dim(a).y != 0) 
            m = get_dim(a).y % get_mem_dim(a).y;
        if(j == (get_grid_dim(a).x-1) && get_dim(a).x%get_mem_dim(a).x != 0) 
            n = get_dim(a).x % get_mem_dim(a).x;
       
        if(m != ld || n != sd){
            for(size_t jj=0; jj<n; jj++){
                for(size_t ii=0; ii<m; ii++)
                    randomize((ad+(jj*ld+ii))); // pointer arithmetic
                if(m != ld) memset(&ad[jj*ld+m], 0, sizeof(T)*(ld-m));
            }
            if(n != sd) 
                memset(&ad[n*ld], 0, sizeof(T)*(sd-n)*ld);
        }else{
            for(size_t ii=0; ii<n*m; ii++)
                randomize((ad+ii)); // pointer arithmetic
        }
    }

    template<typename T>
    void null_i(models::v_model::object& a){
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t m = get_mem_dim(a).y;
        size_t n = get_mem_dim(a).x;
        T* ad = current(a)(i,j);
        memset(ad, 0, m*n*sizeof(T));
    }

    template<typename T>
    void identity_i(models::v_model::object& a){
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t m = get_mem_dim(a).y;
        size_t n = get_mem_dim(a).x;
        T* ad = current(a)(i,j);
        memset(ad, 0, m*n*sizeof(T));
        if((i+1)*m <= j*n) return;
        if(i*m >= (j+1)*n) return;
        for(size_t jj = j*n; jj < (j+1)*n; jj++){
            if(i*m > jj) continue;
            if((i+1)*m <= jj) continue;
            if(jj < get_dim(a).x && jj < get_dim(a).y) 
                ad[jj % m + (jj%n)*m] = 1;
        }
    }

    template<typename T>
    void value_i(models::v_model::object& a){
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t m = get_mem_dim(a).y;
        size_t n = get_mem_dim(a).x;
        T* ad = current(a)(i,j);
        for(size_t ii=0; ii<n*m; ii++) ad[ii] = 0; //a.init_value; // IT DOESN'T EXIST!
    }

}
#endif
