// C - Note about the vectorization, the maximum of iteration is calculated outside the loop, to help at least the intel compiler to vectorize
#include <limits> // for the validation test
namespace ambient {

    #include "ambient/utils/numeric.h" // Blas/Lapack signature
    #include "ambient/utils/ceil.h"

    // C - w, alfa are not nested into dim2 because they are not 2D coordinates
    template <typename T>
    void __a_memcpy(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, maquis::types::p_dense_matrix_impl<T> const& src, T *sd,  dim2 const& dpos, dim2 const& spos, size_t w, double alfa){
        size_t v = get_mem_dim(src).y-spos.x;
        memcpy(&dd[dpos.y*get_mem_dim(dest).y+dpos.x],
               &sd[spos.y*get_mem_dim(src).y+spos.x],
               std::min(v, w)*sizeof(T));
    }

    template <typename T>
    void __a_memscal(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, maquis::types::p_dense_matrix_impl<T> const& src, T *sd,  dim2 const& dpos, dim2 const& spos, size_t w, double alfa){
        size_t v = get_mem_dim(src).y-spos.x;
        for(int z = 0; z < std::min(v, w); z++)
            dd[dpos.y*get_mem_dim(dest).y+dpos.x+z] += sd[spos.y*get_mem_dim(src).y+spos.x + z]*alfa;
    }

    // V = double, S = size_t
    template<typename T>
    void __a_memptf(void (*ptf)(maquis::types::p_dense_matrix_impl<T>& dest, T* dd, 
                                maquis::types::p_dense_matrix_impl<T> const& src, T *sd , 
                                dim2 const& dpos, dim2 const& spos, size_t w, double alfa),
                    maquis::types::p_dense_matrix_impl<T>& dest, dim2 dest_p, 
                    const maquis::types::p_dense_matrix_impl<T>& src, dim2 src_p, 
                    dim2 size, double alfa = 0.0)
    {
        // C - memcopy implementation for ambient - p_dense_matrix representation
        // C - The ouput (dest) must be a pinned p_dense_matrix
        size_t starti, startj, limi, limj;
        size_t di = ctxt.get_block_id().y * get_mem_dim(dest).y;
        size_t dj = ctxt.get_block_id().x * get_mem_dim(dest).x;
    
        assert(get_grid_dim(dest).x*get_work_dim(dest).x - dest_p.x >= size.x);
        assert(get_grid_dim(dest).y*get_work_dim(dest).y - dest_p.y >= size.y);
        assert(get_grid_dim(src).x*get_work_dim(src).x - src_p.x >= size.x);
        assert(get_grid_dim(src).y*get_work_dim(src).y - src_p.y >= size.y);
    
        if(size.x == 0 || size.y == 0) return;
        if((di + get_mem_dim(dest).y <= dest_p.y) || (dj + get_mem_dim(src).x  <= dest_p.x)) return;
        if((di >= dest_p.y + size.y) || (dj >= dest_p.x + size.x)) return;
    // lets find dest-block copy limits
        if(di + get_mem_dim(dest).y > dest_p.y + size.y) limi = (dest_p.y + size.y) % get_mem_dim(dest).y;
        else limi = get_mem_dim(dest).y;
        if(dj + get_mem_dim(dest).x > dest_p.x + size.x) limj = (dest_p.x + size.x) % get_mem_dim(dest).x;
        else limj = get_mem_dim(dest).x;
    // lets find dest-block starting point
        if(di < dest_p.y) starti = dest_p.y % get_mem_dim(dest).y;
        else starti = 0;
        if(dj < dest_p.x) startj = dest_p.x % get_mem_dim(dest).x;
        else startj = 0;
    
        size_t si = di + starti - dest_p.y + src_p.y;
        size_t sii = si % get_mem_dim(src).y;
    // let's find how many blocks do we need for this one
        size_t src_blocks_i = 1;
        int num_src_blocks = limi-starti-get_mem_dim(src).y+sii;
        if(num_src_blocks > 0) src_blocks_i = __a_ceil( num_src_blocks / get_mem_dim(src).y ) + 1;
    // let's exhaust first src block
        T* dd = current(dest)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    
        dim2 dpos,spos;
        for(size_t j = startj; j < limj; j++){
            size_t sj = dj + j - dest_p.x + src_p.x;
            size_t w = limi - starti;
            dpos.x = starti;
            dpos.y = j;
            spos.x = si % get_mem_dim(src).y;
            spos.y = sj % get_mem_dim(src).x;
            for(int k = 0; k < src_blocks_i; k++){
                T* sd = current(src)(si / get_mem_dim(src).y + k, sj / get_mem_dim(src).x);
                ptf(dest,dd,src,sd,dpos,spos,w,alfa);            
                w -= get_mem_dim(src).y-spos.x;
                dpos.x += get_mem_dim(src).y-spos.x;
                spos.x = 0;
            }
        }
    }

    template<typename T>
    void print_pinned_block(T& a){
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
    
        double* ad = current(a)(i, j);
        std::cout << " --(i,j)-- " << i << "," << j << std::endl;
        for(int ii = 0; ii < get_mem_dim(a).y; ii++){
            for(int jj = 0; jj < get_mem_dim(a).x; jj++){
                std::cout << ad[jj*get_mem_dim(a).y + ii]<< " ";
            }     
            std::cout << " " << std::endl;
        }
        std::cout << " --------------------------- " << std::endl;
    }

    template<typename T>
    void gemm_inplace_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |  x  | 0 | 1 | 2 |  =  | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    //
    // partial reduce?..
    /////////////////////////////////////////////////////////////////////////
        int m   = get_mem_dim(a).y;
        int n   = get_mem_dim(b).x;
        int k   = get_mem_dim(b).y;
        int lda = m;
        int ldb = k;
        int ldc = m;
        T alpha(1.0); 
        T beta(1.0);
    // a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
    // current block of matrix a:
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
    // taking (j,i) of b:
        if(get_mem_grid_dim(b).y > j) while(i < get_mem_grid_dim(b).x){
            T* bd = current(b)(j,i); // remote
    // multiplying with column of a:
            for(int z = 0; z < get_mem_grid_dim(a).y; z++){
                T* ad = current(a)(z,j);
                T* cd = reduced<'+'>(a)(z,i); // a(z,j) x b(j,i) => a(z,i)
                gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            }
            i += get_mem_grid_dim(a).y;
        }
    }

    template<typename T>
    void gemm_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |  x  | 0 | 1 | 2 |  =  | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    // | 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
    //  --- --- ---       --- --- ---       --- --- ---
    //
    // partial reduce?..
    /////////////////////////////////////////////////////////////////////////
        int m   = get_mem_dim(a).y;
        int n   = get_mem_dim(b).x;
        int k   = get_mem_dim(b).y;
        int lda = m;
        int ldb = k;
        int ldc = m;
        T alpha(1.0); 
        T beta(1.0);
    // a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
    // current block of matrix a:
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
    // taking (j,i) of b:
        if(get_mem_grid_dim(b).y > j) while(i < get_mem_grid_dim(b).x){
            T* bd = current(b)(j,i); // remote
    // multiplying with column of a:
            for(int z = 0; z < get_mem_grid_dim(a).y; z++){
                T* ad = current(a)(z,j);
                T* cd = updated(c)(z,i); //reduced<'+'>(c)(z,i); // a(z,j) x b(j,i) => c(z,i)
                gemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
            }
            i += get_mem_grid_dim(a).y;
        }
    }

    template<>
    void copy_c(maquis::types::p_dense_matrix_impl<double>& ac, pinned const maquis::types::p_dense_matrix_impl<double>& a){
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        if(i >= get_mem_grid_dim(ac).y || j >= get_mem_grid_dim(ac).x) return;
        double* a_elements  = current(a)(i,j);
        double* ac_elements = updated(ac)(i,j);
        memcpy(ac_elements, a_elements, sizeof(double)*get_mem_dim(a).y*get_mem_dim(a).x);
        //printf("COPYED %d and %d!\n", i, j);
    }

    void variable_free_c(void*& a){ free(a); }

    template<typename T>
    void remove_rows_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& i_mark, const size_t& k){
        // C - Presently I do not copy datas the between num_rows and the lda ....
        // i_mark marks (x position) to remove rows, k number of rows to remove (default 1) 
        size_t numrows = get_dim(a).y;
        size_t numcols = get_dim(a).x;
        __a_memptf(&__a_memcpy<T>, a, dim2(0,i_mark), a, dim2(0,k+i_mark), dim2(numcols,numrows-k-i_mark));
    }

    template<typename T>
    void remove_cols_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& j_mark, const size_t& k){
        // C - Presently I do not copy datas the between num_cols and the sda ....
        // j_mark marks (x position) to remove cols, k number of columns to remove (default 1) 
        size_t numrows = get_dim(a).y;
        size_t numcols = get_dim(a).x;
        __a_memptf(&__a_memcpy<T>,a, dim2(j_mark,0), a, dim2(k+j_mark,0), dim2(numcols-k-j_mark,numrows));
    }

    template<typename T>
    void resize_c(maquis::types::p_dense_matrix_impl<T>& a, const size_t& rows, const size_t& cols){
        __a_memptf(&__a_memcpy<T>, a, dim2(0,0), a, dim2(0,0), dim2(rows,cols));
    }

    template<typename T>
    void sqrt_diagonal_c(pinned maquis::types::p_dense_matrix_impl<T>& a){
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = get_mem_dim(a).y;
        for(int i=0; i < size; i++)
            ad[i] = sqrt(ad[i]);
    }

    template<typename T>
    void exp_diagonal_c(pinned maquis::types::p_dense_matrix_impl<T>& a){
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = get_mem_dim(a).y;
        for(int i=0; i < size; i++)
            ad[i] = exp(ad[i]);
    }

    template<typename T>
    void copy_after_c(pinned maquis::types::p_dense_matrix_impl<T>& ac, const size_t& pos, const maquis::types::p_dense_matrix_impl<T>& a){
        __a_memptf(&__a_memcpy<double>, ac, dim2(0,pos), a, dim2(0,0), dim2(1,get_dim(a).y));
    }

    void copy_after_std_c(std::vector<double>*& ac, const size_t& pos, pinned const maquis::types::p_dense_matrix_impl<double>& a){ // C - bug if execution independant
        double* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = std::min(ac->size(),pos+get_mem_dim(a).y);
        for(int i=0; (pos+i) < size; i++){
            (*ac)[pos+i] = ad[i];
        }
    }

    void push_back_sqr_gt_c(std::vector<double>*& ac, pinned const maquis::types::p_dense_matrix_impl<double>& a, const double& prec){
    double* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    for(int i=0; i < get_mem_dim(a).y; i++)
        if(ad[i] > prec) ac->push_back(ad[i]*ad[i]);
    }

    template<typename T>
    void cast_to_dense_c(std::vector<T>*& ac, pinned const maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n){
        // note as we cast to dense matrix, l_kernel on 1 proc, thus half work done
        int offset,size_y(get_mem_dim(a).y),size_x(get_mem_dim(a).x);
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        int yi = get_mem_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
        int xj = get_mem_dim(a).x*j; 
        T* ad = current(a)(i,j);
       
        //operator ?, case 1 matrix is lower than one work group, case 2 several work groups or fit in x direction
        if(j+1 == get_mem_grid_dim(a).x)
            size_x = (get_mem_dim(a).x > n) ? n : (n - (get_mem_grid_dim(a).x-1)*get_mem_dim(a).x);
    
        for(int ii=0; ii < size_x; ++ii){
            offset = yi + (xj+ii)*m;
            if(i+1 == get_mem_grid_dim(a).y)
                size_y = (get_mem_dim(a).y > m) ? m : (m - (get_mem_grid_dim(a).y-1)*get_mem_dim(a).y); // y direction
    
            memcpy((void*)&(*ac)[offset],(void*)&ad[ii*get_mem_dim(a).x], size_y*sizeof(T));  
        }
    }

    template <typename T>
    void cast_to_p_dense_c(const std::vector<T>*& ac, pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& m, const size_t& n, const size_t& lda){
        int offset,size_y(get_mem_dim(a).y),size_x(get_mem_dim(a).x);
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        int yi = get_mem_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
        int xj = get_mem_dim(a).x*j; 
        T* ad = updated(a)(i,j);
    
        //operator ?, case 1 matrix is lower than one work group, case 2 several work groups or fit in x direction
        if(j+1 == get_mem_grid_dim(a).x)
            size_x = (get_mem_dim(a).x > n) ? n : (n - (get_mem_grid_dim(a).x-1)*get_mem_dim(a).x);
    
        for(int ii=0; ii < size_x; ++ii){
            offset = yi + (xj+ii)*lda; //lda because possible resize;
            if(i+1 == get_mem_grid_dim(a).y)
               size_y = (get_mem_dim(a).y > m) ? m : (m - (get_mem_grid_dim(a).y-1)*get_mem_dim(a).y);
            memcpy((void*)&ad[ii*get_mem_dim(a).x],(void*)&(*ac)[offset], size_y*sizeof(T)); // y direction 
        }
    }

    template <typename T>
    void reshape_l2r_c(const maquis::types::p_dense_matrix_impl<T>& left, pinned maquis::types::p_dense_matrix_impl<T>& right,
                       const size_t& left_offset, const size_t& right_offset, 
                       const size_t& sdim, const size_t& ldim, const size_t& rdim)
    {
        //printf("reshape_l2r_c\n");
        for(size_t ss = 0; ss < sdim; ++ss)
            __a_memptf(&__a_memcpy<T>, right, dim2(ss*rdim + right_offset,0), 
                       left,  dim2(0, ss*ldim + left_offset), 
                       dim2( rdim, ldim ));
    }

    template <typename T>
    void reshape_r2l_c(pinned maquis::types::p_dense_matrix_impl<T>& left, const maquis::types::p_dense_matrix_impl<T>& right,
                       const size_t& left_offset, const size_t& right_offset, 
                       const size_t& sdim, const size_t& ldim, const size_t& rdim)
    {
        for(size_t ss = 0; ss < sdim; ++ss)
            __a_memptf(&__a_memcpy<T>, left,  dim2(0, ss*ldim + left_offset), 
                       right, dim2(ss*rdim + right_offset,0), 
                       dim2( rdim, ldim ));
    }

    template<typename T>
    void __a_add_scaled(T& dest, dim2 dest_p, const T& src, dim2 src_p, typename T::value_type alfa, dim2 size){
        __a_memptf(&__a_memscal<double>, dest, dest_p, src, src_p, size, alfa);
    }

    template <typename T>
    void rb_tensor_mpo_c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                         const size_t& out_offset, const size_t& in_offset, 
                         const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
    {
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                T* alfad = current(alfa)(ss1/get_mem_dim(alfa).y, ss2/get_mem_dim(alfa).x);
                T  alfa_t = alfad[ss1%get_mem_dim(alfa).y + get_mem_dim(alfa).y*(ss2%get_mem_dim(alfa).x)];
                __a_add_scaled(out, dim2(out_offset + ss2*rdim, 0),
                               in,  dim2(in_offset + ss1*rdim, 0),
                               alfa_t, dim2(rdim, ldim));
            }
    }

    template <typename T>
    void lb_tensor_mpo_c(pinned maquis::types::p_dense_matrix_impl<T>& out, const maquis::types::p_dense_matrix_impl<T>& in, const maquis::types::p_dense_matrix_impl<T>& alfa,
                         const size_t& out_offset, const size_t& in_offset, 
                         const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
    {
        for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
            for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
                T* alfad = current(alfa)(ss1/get_mem_dim(alfa).y, ss2/get_mem_dim(alfa).x);
                T  alfa_t = alfad[ss1%get_mem_dim(alfa).y + get_mem_dim(alfa).y*(ss2%get_mem_dim(alfa).x)];
                __a_add_scaled(out, dim2(0, out_offset + ss2*ldim),
                               in,  dim2(0, in_offset + ss1*ldim),
                               alfa_t, dim2(rdim, ldim));
            }
    }

    template<typename T>
    void scalar_norm_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, T*& norm){
        T summ = 0;
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = current(a)(i,j);
        size_t size = get_mem_dim(a).x*get_mem_dim(a).y;
        for(size_t ii=0; ii < size; ii++)
            summ += ad[ii]*ad[ii];
    
        *norm += summ;
    }

    template<typename T>
    void scalar_overlap_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, T*& overlap){
        T summ = 0;
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = current(a)(i,j);
        T* bd = current(b)(i,j);
        size_t size = get_mem_dim(a).x*get_mem_dim(a).y;
        for(size_t ii=0; ii < size; ii++)
            summ += ad[ii]*bd[ii];
        *overlap += summ;
    }

    template<typename T>
    void add_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = current(a)(i, j);
        T* bd = current(b)(i, j);
        T* ar = updated(a)(i, j);
        size_t size = get_mem_dim(a).x*get_mem_dim(a).y;
        for(int k = 0; k < size; k++)
            ar[k] = ad[k] + bd[k];
        //printf("ADDED %d and %d! (%lu)\n", i, j, size);
    }

    template<typename T>
    void sub_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b){
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        T* bd = current(b)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = get_mem_dim(a).x*get_mem_dim(a).y;
        for(int i=0; i < size; i++)
            ad[i] -= bd[i];
    }

    template<typename T, typename T2>
    void scale_c(pinned maquis::types::p_dense_matrix_impl<T>& m, const T2& t){
        T* md   = current(m)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        size_t size = get_mem_dim(m).x*get_mem_dim(m).y;
        for(int i=0; i < size; i++)
            md[i] *= t;
    }

    template<typename T>
    void gemm_diagonal_lhs_c(const maquis::types::p_dense_matrix_impl<T>& a_diag, pinned const maquis::types::p_dense_matrix_impl<T>& b, maquis::types::p_dense_matrix_impl<T>& c){
        int j = ctxt.get_block_id().y*get_mem_dim(b).y;
        int size = get_mem_dim(b).x;
        int lda  = get_mem_dim(b).y;
        int ONE  = 1;
        T* bd = current(b)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        T* cd = current(c)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    
        memset(cd, 0, get_mem_dim(c).x*get_mem_dim(c).y*sizeof(T));
        for(int jj = 0 ; jj < get_mem_dim(b).y ; jj++){
             T* alpha = current(a_diag)((j+jj)/get_mem_dim(a_diag).y,0);
    	 axpy(&size, &alpha[(j+jj)%get_mem_dim(a_diag).y], &bd[jj], &lda, &cd[jj], &lda);
        }
    }

    template<typename T>
    void gemm_diagonal_rhs_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b_diag, maquis::types::p_dense_matrix_impl<T>& c){
        //printf("gemm_diagonal_rhs_c\n");
        int j = ctxt.get_block_id().x*get_mem_dim(a).x;
        int size = get_mem_dim(a).y;
        int ONE = 1;
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        T* cd = current(c)(ctxt.get_block_id().y, ctxt.get_block_id().x);
    
        memset(cd, 0, get_mem_dim(c).x*get_mem_dim(c).y*sizeof(T));
        for(int jj = 0 ; jj < get_mem_dim(a).x ; jj++){
    	 T* alpha = current(b_diag)((j+jj)/get_mem_dim(b_diag).y,0);
    	 axpy(&size, &alpha[(j+jj)%get_mem_dim(b_diag).y], &ad[jj*get_mem_dim(a).y], &ONE, &cd[jj*get_mem_dim(c).y], &ONE);
        }
    }

    template<typename T>
    void trace_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, T*& trace){ // originated from identity_i
      // TO DO : O think this trace kernel only works if the a fits into one workgroup
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        size_t m = get_mem_dim(a).y;
        size_t n = get_mem_dim(a).x;
        T* ad = current(a)(i,j);
    
        if((i+1)*m <= j*n) return;
        if(i*m >= (j+1)*n) return;
        for(size_t jj = j*n; jj < std::min((j+1)*n,(size_t)get_dim(a).x); jj++){
            if(i*m > jj) continue;
            if((i+1)*m <= jj) continue;
           *trace += ad[jj % m + (jj%n)*m];
        }
    }

    template<typename T>
    void transpose_c(pinned maquis::types::p_dense_matrix_impl<T>& transposed, const maquis::types::p_dense_matrix_impl<T>& original){
        size_t i = ctxt.get_block_id().y;
        size_t j = ctxt.get_block_id().x;
        T* td = current(transposed)(i,j);
        T* od = current(original)(j,i);
    
        for(size_t i = 0; i < get_mem_dim(original).y; ++i){
            for(size_t j=0; j < get_mem_dim(original).x; ++j){
                td[j+i*get_mem_dim(transposed).y] = od[i+j*get_mem_dim(original).y];
            }
        }
    }

    template<typename T>
    void nullcut_c(pinned maquis::types::p_dense_matrix_impl<T>& a, const size_t& num_rows, const size_t& num_cols){
        size_t i = ctxt.get_block_id().y*get_mem_dim(a).y; 
        size_t j = ctxt.get_block_id().x*get_mem_dim(a).x; 
        if((i+get_mem_dim(a).y <= num_rows) && (j+get_mem_dim(a).x <= num_cols)) return;
    
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        for(size_t jj = 0; jj < get_mem_dim(a).x; jj++){
            if(j+jj < num_cols && (i+get_mem_dim(a).y > num_rows)) memset(&ad[jj*get_mem_dim(a).y+num_rows%get_mem_dim(a).y], 0, (get_mem_dim(a).y-num_rows%get_mem_dim(a).y)*sizeof(T));
            else if(j+jj >= num_cols) memset(&ad[jj*get_mem_dim(a).y], 0, get_mem_dim(a).y*sizeof(T));
        }
    }

    template <typename T> 
    void validation_c(pinned const maquis::types::p_dense_matrix_impl<T>& a, const maquis::types::p_dense_matrix_impl<T>& b, int*& bl){ // see paper for Reference Dongara 
        int i = ctxt.get_block_id().y;
        int j = ctxt.get_block_id().x;
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x); 
        T* bd = current(b)(ctxt.get_block_id().y, ctxt.get_block_id().x); 
        double res(0.0); 
        double epsilon = std::numeric_limits<double>::epsilon();
        int position_x(0),position_y(0),position_xy(0); 
    
        for(int ii=0; ii < get_mem_dim(a).x; ++ii){ // the std::resize + cast dense to p_dense can add ghost elements between lda + num_rows
            for(int jj=0; jj < get_mem_dim(a).y; ++jj){
                position_x = j*get_mem_dim(a).x+jj;
                position_y = i*get_mem_dim(a).x+ii;
                if(position_x < get_dim(a).x && position_y < get_dim(a).y){
                    position_xy = ii*get_mem_dim(a).y+jj;
                    res = (norm(ad[position_xy])-norm(bd[position_xy]))/fabs(epsilon*norm(bd[position_xy])); // to do : rotation pb  with complex to change
                    if(res > 256){ // 16 is recommended by Dongara, 256 because lapack gives != runs after runs
                        std::cout <<   ctxt.get_block_id().y << " " <<  ctxt.get_block_id().x << " " << res << " " << ad[i] << " " << bd[i] << std::endl; // C - cout because double or complex
                    *bl = 0; // test failed return 0 (bool false)
                    }
                }
            }
        } 
    }

    // {{{ MKL LAPACK kernels

    template<typename T>
    void svd_c(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n, maquis::types::p_dense_matrix_impl<T>& u, maquis::types::p_dense_matrix_impl<T>& vt, maquis::types::p_dense_matrix_impl<double>& s){
    /* Locals */
        int lda = get_grid_dim(a).y*get_work_dim(a).y;
        int ldu = get_grid_dim(u).y*get_work_dim(u).y;
        int ldvt = get_grid_dim(vt).y*get_work_dim(vt).y;
        int info, lwork;
        T wkopt;
        T* work;
        double* rwork = new double[5*std::min(lda,ldu)]; // C - useless for double but need for complex 
        T* ad = (T*)solidify(current(a));
        T* ud = (T*)solidify(current(u));
        T* vtd = (T*)solidify(current(vt));
        double* sd = (T*)solidify(current(s));
    /* Query and allocate the optimal workspace */
        lwork = -1; // C - Alex, netlib said -1 for the best workspace
        gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, &wkopt, &lwork, rwork, &info );
        lwork = OptimalSize(wkopt);
        work = (T*)malloc( lwork*sizeof(T) );
    /* Compute SVD */
        gesvd( "S", "S", &m, &n, ad, &lda, sd, ud, &ldu, vtd, &ldvt, work, &lwork, rwork, &info );
    /* Check for convergence */
        if( info > 0 ) {
            printf( "The algorithm computing SVD failed to converge.\n" );
            exit( 1 );
        }
        disperse(ud, current(u));
        disperse(vtd, current(vt));
        disperse(sd, current(s));
        free(work);
    }

    void syev_c(const maquis::types::p_dense_matrix_impl<double>& a, int& m, maquis::types::p_dense_matrix_impl<double>& w){
         int lda = get_grid_dim(a).y*get_work_dim(a).y;
         int info, lwork = -1;
    
         double wkopt;
         double* work;
         double* ad = (double*)solidify(current(a));
         double* wd = (double*)solidify(current(w));
         
         dsyev_("V","U",&m,ad,&lda,wd,&wkopt,&lwork,&info);
         lwork = (int)wkopt;
         work = (double*)malloc( lwork*sizeof(double) );
         dsyev_("V","U",&m,ad,&lda,wd,work,&lwork,&info);
    
         if( info > 0 ) {
             printf( "The algorithm computing SYEV failed to converge.\n" );
             exit( 1 );
         }
    
         disperse(ad, current(a));
         disperse(wd, current(w));
         free(work); 
    }

    void heev_c(const maquis::types::p_dense_matrix_impl<double>& a, int& m, maquis::types::p_dense_matrix_impl<double>& w){
         int lda = get_grid_dim(a).y*get_work_dim(a).y;
         int info, lwork = -1;
    
         double wkopt;
         double* work;
         double* ad = (double*)solidify(current(a));
         double* wd = (double*)solidify(current(w));
    
         dsyev_("V","U",&m,ad,&lda,wd,&wkopt,&lwork,&info);
         lwork = (int)wkopt;
         work = (double*)malloc( lwork*sizeof(double) );
         dsyev_("V","U",&m,ad,&lda,wd,work,&lwork,&info);
    
         if( info > 0 ) {
             printf( "The algorithm computing SYEV failed to converge.\n" );
             exit( 1 );
         }
        
         // First we reverse the eigenvalues, to be in agreement with the serial version ! 
         // The matrix is solidified, so we do not care on the workgroup representation
         double tempdbl;
         for (int i=0; i< static_cast<int>(m/2); i++){ 
             tempdbl = wd[i];
             wd[i] = wd[m-i-1];
             wd[m-i-1] = tempdbl;
         } 
         // Second we reverse the eigenvectors
         double* tempcol = new double[lda]; 
         for (int i=0; i< static_cast<int>(m/2); ++i){ 
             memmove((void*)tempcol,(void*)&ad[i*lda],lda*sizeof(double));
             memmove((void*)&ad[i*lda],(void*)&ad[(m-1-i)*lda],lda*sizeof(double));
             memmove((void*)&ad[(m-1-i)*lda],(void*)tempcol,lda*sizeof(double));
         }
         delete[] tempcol; 
     
         disperse(ad, current(a));
         disperse(wd, current(w));
         free(work);
    }

    // }}}

    template<typename T>
    void print_c(const maquis::types::p_dense_matrix_impl<T>& a, int& m, int& n){
        int lda = get_grid_dim(a).y*get_work_dim(a).y;
        T* pa = (T*)solidify(current(a));
        
        for(int i=0; i < m; ++i){
            for(int j=0; j < n; ++j){
                std::cout << *(pa+i+j*lda) << " ";                    
            }
            printf("\n");
        }
    }

    template<typename T>
    void initv_c(pinned maquis::types::p_dense_matrix_impl<T>&a, T const&v){ // for devs
        T* ad = current(a)(ctxt.get_block_id().y, ctxt.get_block_id().x);
        for(int i=0; i < get_mem_dim(a).x*get_mem_dim(a).y; ++i)
            ad[i] = v;
    }

}
