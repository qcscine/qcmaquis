// C - Note about the vectorization, the maximum of iteration  is calculated outside the for loop, to help the intel compiler to vectorize

#include "numeric.h" // Blas/Lapack signature

#define NODE_COUNT 1

// C - w, alfa are not nested into dim2 because they are not 2D coordinates
template <typename T, typename V>
void __a_memcpy(T& dest, V* dd, T const& src, V *sd,  dim2 const& dpos, dim2 const& spos, std::size_t w, double alfa){
    std::size_t v = get_mem_t_dim(src).y-spos.x;
    memcpy(&dd[dpos.y*get_mem_t_dim(dest).y+dpos.x],
           &sd[spos.y*get_mem_t_dim(src).y+spos.x],
           std::min(v, w)*sizeof(double));
}

template <typename T, typename V>
void __a_memscal(T& dest, V* dd, T const& src, V *sd,  dim2 const& dpos, dim2 const& spos, std::size_t w, double alfa){
    std::size_t v = get_mem_t_dim(src).y-spos.x;
    for(int z = 0; z < std::min(v, w); z++)
        dd[dpos.y*get_mem_t_dim(dest).y+dpos.x+z] += sd[spos.y*get_mem_t_dim(src).y+spos.x + z]*alfa;
}
// V = double, S = size_t
template<typename T, typename V, typename S>
void __a_memptf(void (*ptf)(T& dest, V* dd, T const& src, V *sd , dim2 const& dpos, dim2 const& spos, S w, V alfa),
                T& dest, dim2 dest_p, const T& src, dim2 src_p, dim2 size, V alfa = 0.0)
{
    // C - memcopy implementation for ambient - p_dense_matrix representation
    // C - The ouput (dest) must be a pinned p_dense_matrix
    size_t starti, startj, limi, limj;
    size_t di = get_block_id(dest).y * get_mem_t_dim(dest).y;
    size_t dj = get_block_id(dest).x * get_mem_t_dim(dest).x;

    assert(get_grid_dim(dest).x*get_mem_t_dim(dest).x - dest_p.x >= size.x);
    assert(get_grid_dim(dest).y*get_mem_t_dim(dest).y - dest_p.y >= size.y);
    assert(get_grid_dim(src).x*get_mem_t_dim(src).x - src_p.x >= size.x);
    assert(get_grid_dim(src).y*get_mem_t_dim(src).y - src_p.y >= size.y);

    if(size.x == 0 || size.y == 0) return;
    if((di + get_mem_t_dim(dest).y <= dest_p.y) || (dj + get_mem_t_dim(src).x  <= dest_p.x)) return;
    if((di >= dest_p.y + size.y) || (dj >= dest_p.x + size.x)) return;
// lets find dest-block copy limits
    if(di + get_mem_t_dim(dest).y > dest_p.y + size.y) limi = (dest_p.y + size.y) % get_mem_t_dim(dest).y;
    else limi = get_mem_t_dim(dest).y;
    if(dj + get_mem_t_dim(dest).x > dest_p.x + size.x) limj = (dest_p.x + size.x) % get_mem_t_dim(dest).x;
    else limj = get_mem_t_dim(dest).x;
// lets find dest-block starting point
    if(di < dest_p.y) starti = dest_p.y % get_mem_t_dim(dest).y;
    else starti = 0;
    if(dj < dest_p.x) startj = dest_p.x % get_mem_t_dim(dest).x;
    else startj = 0;

    size_t si = di + starti - dest_p.y + src_p.y;
    size_t sii = si % get_mem_t_dim(src).y;
// let's find how many blocks do we need for this one
    size_t src_blocks_i = 1;
    int num_src_blocks = limi-starti-get_mem_t_dim(src).y+sii;
    if(num_src_blocks > 0) src_blocks_i = __a_ceil( num_src_blocks / get_mem_t_dim(src).y ) + 1;
// let's exhaust first src block
    typename T::value_type* dd = current(dest)(get_block_id(dest).y, get_block_id(dest).x);

    dim2 dpos,spos;
    for(size_t j = startj; j < limj; j++){
        size_t sj = dj + j - dest_p.x + src_p.x;
        size_t w = limi - starti;
        dpos.x = starti;
        dpos.y = j;
        spos.x = si % get_mem_t_dim(src).y;
        spos.y = sj % get_mem_t_dim(src).x;
        for(int k = 0; k < src_blocks_i; k++){
            typename T::value_type* sd = current(src)(si / get_mem_t_dim(src).y + k, sj / get_mem_t_dim(src).x);
            ptf(dest,dd,src,sd,dpos,spos,w,alfa);            
            w -= get_mem_t_dim(src).y-spos.x;
            dpos.x += get_mem_t_dim(src).y-spos.x;
            spos.x = 0;
        }
    }
}

template<typename T>
void print_pinned_block(T& a)
{
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;

    double* ad = current(a)(i, j);
    std::cout << " --(i,j)-- " << i << "," << j << std::endl;
    for(int ii = 0; ii < get_mem_t_dim(a).y; ii++){
        for(int jj = 0; jj < get_mem_t_dim(a).x; jj++){
            std::cout << ad[jj*get_mem_t_dim(a).y + ii]<< " ";
        }     
        std::cout << " " << std::endl;
    }
    std::cout << " --------------------------- " << std::endl;
}

void gemm_inplace_c(pinned p_dense_matrix<double>& a, const p_dense_matrix<double>& b)
{
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
    int m   = get_mem_t_dim(a).y;
    int n   = get_mem_t_dim(b).x;
    int k   = get_mem_t_dim(b).y;
    int lda = m;
    int ldb = k;
    int ldc = m;
    double alpha = 1.0; 
    double beta  = 1.0;
// a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
// current block of matrix a:
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
// taking (j,i) of b:
    if(get_grid_dim(b).y > j) while(i < get_grid_dim(b).x){
        double* bd = current(b)(j,i); // remote
// multiplying with column of a:
        for(int z = 0; z < get_grid_dim(a).y; z++){
            double* ad = current(a)(z,j);
            double* cd = reduced<'+'>(a)(z,i); // a(z,j) x b(j,i) => a(z,i)
            dgemm_("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
        }
        i += get_grid_dim(a).y;
    }
}

void gemm_c(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
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
    int m   = get_mem_t_dim(a).y;
    int n   = get_mem_t_dim(b).x;
    int k   = get_mem_t_dim(b).y;
    int lda = m;
    int ldb = k;
    int ldc = m;
    double alpha = 1.0; 
    double beta  = 1.0;
// a(i,j) => a(z,j) x b(j,i)  where z : [0,m)
// current block of matrix a:
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
// taking (j,i) of b:
    if(get_grid_dim(b).y > j) while(i < get_grid_dim(b).x){
        double* bd = current(b)(j,i); // remote
// multiplying with column of a:
        for(int z = 0; z < get_grid_dim(a).y; z++){
            double* ad = current(a)(z,j);
            double* cd = reduced<'+'>(c)(z,i); // a(z,j) x b(j,i) => c(z,i)
            dgemm_("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
        }
        i += get_grid_dim(a).y;
    }
}

void copy_c(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    if(i >= get_grid_dim(ac).y || j >= get_grid_dim(ac).x) return;
    double* a_elements  = current(a)(i,j);
    double* ac_elements = current(ac)(i,j);
    memcpy(ac_elements, a_elements, sizeof(double)*get_mem_t_dim(a).y*get_mem_t_dim(a).x);
}

void touch_c(const p_dense_matrix<double>& a){ }

void variable_free_c(void*& a){ free(a); }

void remove_rows_c(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    // C - Presently I do not copy datas the between num_rows and the lda ....
    // i_mark marks (x position) to remove rows, k number of rows to remove (default 1) 
    std::size_t numrows = get_dim(a).y;
    std::size_t numcols = get_dim(a).x;
    __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>, a, dim2(0,i_mark), a, dim2(0,k+i_mark), dim2(numcols,numrows-k-i_mark));
}

void remove_cols_c(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{
    // C - Presently I do not copy datas the between num_cols and the sda ....
    // j_mark marks (x position) to remove cols, k number of columns to remove (default 1) 
    std::size_t numrows = get_dim(a).y;
    std::size_t numcols = get_dim(a).x;
    __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>,a, dim2(j_mark,0), a, dim2(k+j_mark,0), dim2(numcols-k-j_mark,numrows));
}

void resize_c(p_dense_matrix<double>& a, const size_t& rows, const size_t& cols)
{
    __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>, a, dim2(0,0), a, dim2(0,0), dim2(rows,cols));
}

void sqrt_diagonal_c(pinned p_dense_matrix<double>& a)
{
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    std::size_t size = get_mem_t_dim(a).y;
    for(int i=0; i < size; i++)
        ad[i] = sqrt(ad[i]);
}

void exp_diagonal_c(pinned p_dense_matrix<double>& a)
{
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    std::size_t size = get_mem_t_dim(a).y;
    for(int i=0; i < size; i++)
        ad[i] = exp(ad[i]);
}

void copy_after_c(pinned p_dense_matrix<double>& ac, const size_t& pos, const p_dense_matrix<double>& a)
{
    __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>, ac, dim2(0,pos), a, dim2(0,0), dim2(1,get_dim(a).y));
}

void copy_after_std_c(std::vector<double>*& ac, const size_t& pos, pinned const p_dense_matrix<double>& a)
{ // C - bug if execution independant
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    std::size_t size = std::min(ac->size(),pos+get_mem_t_dim(a).y);
    for(int i=0; (pos+i) < size; i++){
        (*ac)[pos+i] = ad[i];
    }
}

void push_back_sqr_gt_c(std::vector<double>*& ac, pinned const p_dense_matrix<double>& a, const double& prec)
{
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    for(int i=0; i < get_mem_t_dim(a).y; i++)
        if(ad[i] > prec) ac->push_back(ad[i]*ad[i]);
}

void cast_to_dense_c(std::vector<double>*& ac, pinned const p_dense_matrix<double>& a, const size_t& m, const size_t& n)
{ // note as we cast to dense matrix, l_kernel on 1 proc, thus half work done
    int offset,size_y(get_mem_t_dim(a).y),size_x(get_mem_t_dim(a).x);
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    int yi = get_mem_t_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
    int xj = get_mem_t_dim(a).x*j; 
    double* ad = current(a)(i,j);
   
    //operator ?, case 1 matrix is lower than one work group, case 2 several work groups or fit in x direction
    if(j+1 == get_grid_dim(a).x)
        size_x = (get_mem_t_dim(a).x > n) ? n : (n - (get_grid_dim(a).x-1)*get_mem_t_dim(a).x);

    for(int ii=0; ii < size_x; ++ii){
        offset = yi + (xj+ii)*m;
        if(i+1 == get_grid_dim(a).y)
            size_y = (get_mem_t_dim(a).y > m) ? m : (m - (get_grid_dim(a).y-1)*get_mem_t_dim(a).y); // y direction

        memcpy((void*)&(*ac)[offset],(void*)&ad[ii*get_mem_t_dim(a).x], size_y*sizeof(double));  
    }
}

void cast_to_p_dense_c(const std::vector<double>*& ac, pinned p_dense_matrix<double>& a, const size_t& m, const size_t& n, const size_t& lda)
{
    int offset,size_y(get_mem_t_dim(a).y),size_x(get_mem_t_dim(a).x);
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    double* ad = current(a)(i,j);
    int yi = get_mem_t_dim(a).y*i; // conversion cartersia coordinates dense / p_dense
    int xj = get_mem_t_dim(a).x*j; 

    //operator ?, case 1 matrix is lower than one work group, case 2 several work groups or fit in x direction
    if(j+1 == get_grid_dim(a).x)
        size_x = (get_mem_t_dim(a).x > n) ? n : (n - (get_grid_dim(a).x-1)*get_mem_t_dim(a).x);

    for(int ii=0; ii < size_x; ++ii){
        offset = yi + (xj+ii)*lda; //lda because possible resize;
        if(i+1 == get_grid_dim(a).y)
           size_y = (get_mem_t_dim(a).y > m) ? m : (m - (get_grid_dim(a).y-1)*get_mem_t_dim(a).y);
        memcpy((void*)&ad[ii*get_mem_t_dim(a).x],(void*)&(*ac)[offset], size_y*sizeof(double)); // y direction 
    }
}

void reshape_l2r_c(const p_dense_matrix<double>& left, pinned p_dense_matrix<double>& right,
                   const size_t& left_offset, const size_t& right_offset, 
                   const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    //printf("reshape_l2r_c\n");
    for(size_t ss = 0; ss < sdim; ++ss)
        __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>, right, dim2(ss*rdim + right_offset,0), 
                   left,  dim2(0, ss*ldim + left_offset), 
                   dim2( rdim, ldim ));
}

void reshape_r2l_c(pinned p_dense_matrix<double>& left, const p_dense_matrix<double>& right,
                   const size_t& left_offset, const size_t& right_offset, 
                   const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    for(size_t ss = 0; ss < sdim; ++ss)
        __a_memptf(&__a_memcpy<p_dense_matrix<double>, double>, left,  dim2(0, ss*ldim + left_offset), 
                   right, dim2(ss*rdim + right_offset,0), 
                   dim2( rdim, ldim ));
}

template<typename T>
void __a_add_scaled(T& dest, dim2 dest_p, const T& src, dim2 src_p, typename T::value_type alfa, dim2 size)
{
    __a_memptf(&__a_memscal<p_dense_matrix<double>, double>, dest, dest_p, src, src_p, size, alfa);
}

template <typename T>
void rb_tensor_mpo_c(pinned p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                     const size_t& out_offset, const size_t& in_offset, 
                     const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
{
    for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
            T* alfad = current(alfa)(ss1/get_mem_t_dim(alfa).y, ss2/get_mem_t_dim(alfa).x);
            T  alfa_t = alfad[ss1%get_mem_t_dim(alfa).y + get_mem_t_dim(alfa).y*(ss2%get_mem_t_dim(alfa).x)];
            __a_add_scaled(out, dim2(out_offset + ss2*rdim, 0),
                           in,  dim2(in_offset + ss1*rdim, 0),
                           alfa_t, dim2(rdim, ldim));
        }
}

template <typename T>
void lb_tensor_mpo_c(pinned p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                     const size_t& out_offset, const size_t& in_offset, 
                     const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
{
    for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
            T* alfad = current(alfa)(ss1/get_mem_t_dim(alfa).y, ss2/get_mem_t_dim(alfa).x);
            T  alfa_t = alfad[ss1%get_mem_t_dim(alfa).y + get_mem_t_dim(alfa).y*(ss2%get_mem_t_dim(alfa).x)];
            __a_add_scaled(out, dim2(0, out_offset + ss2*ldim),
                           in,  dim2(0, in_offset + ss1*ldim),
                           alfa_t, dim2(rdim, ldim));
        }
}

void scalar_norm_c(pinned const p_dense_matrix<double>& a, double*& norm)
{
    double summ = 0;
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    double* ad = current(a)(i,j);
    std::size_t size = get_mem_t_dim(a).x*get_mem_t_dim(a).y;
    for(size_t ii=0; ii < size; ii++)
        summ += ad[ii]*ad[ii];

    *norm += summ;
}

void scalar_overlap_c(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, double*& overlap)
{
    double summ = 0;
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    double* ad = current(a)(i,j);
    double* bd = current(b)(i,j);
    std::size_t size = get_mem_t_dim(a).x*get_mem_t_dim(a).y;
    for(size_t ii=0; ii < size; ii++)
        summ += ad[ii]*bd[ii];
    *overlap += summ;
}

template<typename T>
void add_c(pinned p_dense_matrix<T>& a, const p_dense_matrix<T>& b)
{
    T* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    T* bd = current(b)(get_block_id(a).y, get_block_id(a).x);
    std::size_t size = get_mem_t_dim(a).x*get_mem_t_dim(a).y;
    for(int i=0; i < size; i++)
        ad[i]+=bd[i];
}

template<typename T>
void sub_c(pinned p_dense_matrix<T>& a, const p_dense_matrix<T>& b)
{
    T* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    T* bd = current(b)(get_block_id(a).y, get_block_id(a).x);
    std::size_t size = get_mem_t_dim(a).x*get_mem_t_dim(a).y;
    for(int i=0; i < size; i++)
        ad[i] -= bd[i];
}

template<typename T, typename T2>
void scale_c(pinned p_dense_matrix<T>& m, const T2& t)
{
    T* md   = current(m)(get_block_id(m).y, get_block_id(m).x);
    std::size_t size = get_mem_t_dim(m).x*get_mem_t_dim(m).y;
    for(int i=0; i < size; i++)
        md[i] *= t;
}

void gemm_diagonal_lhs_c(const p_dense_matrix<double>& a_diag, pinned const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    int j = get_block_id(b).y*get_mem_t_dim(b).y;
    int size = get_mem_t_dim(b).x;
    int lda  = get_mem_t_dim(b).y;
    int ONE  = 1;
    double* bd = current(b)(get_block_id(b).y, get_block_id(b).x);
    double* cd = current(c)(get_block_id(b).y, get_block_id(b).x);

    memset(cd, 0, get_mem_t_dim(c).x*get_mem_t_dim(c).y*sizeof(double));
    for(int jj = 0 ; jj < get_mem_t_dim(b).y ; jj++){
         double* alpha = current(a_diag)((j+jj)/get_mem_t_dim(a_diag).y,0);
	 daxpy_(&size, &alpha[(j+jj)%get_mem_t_dim(a_diag).y], &bd[jj], &lda, &cd[jj], &lda);
    }
}

void gemm_diagonal_rhs_c(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b_diag, p_dense_matrix<double>& c)
{
    //printf("gemm_diagonal_rhs_c\n");
    int j = get_block_id(a).x*get_mem_t_dim(a).x;
    int size = get_mem_t_dim(a).y;
    int ONE = 1;
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    double* cd = current(c)(get_block_id(a).y, get_block_id(a).x);

    memset(cd, 0, get_mem_t_dim(c).x*get_mem_t_dim(c).y*sizeof(double));
    for(int jj = 0 ; jj < get_mem_t_dim(a).x ; jj++){
	 double* alpha = current(b_diag)((j+jj)/get_mem_t_dim(b_diag).y,0);
	 daxpy_(&size, &alpha[(j+jj)%get_mem_t_dim(b_diag).y], &ad[jj*get_mem_t_dim(a).y], &ONE, &cd[jj*get_mem_t_dim(c).y], &ONE);
    }
}

void trace_c(pinned const p_dense_matrix<double>& a, double*& trace)
{ // originated from identity_i
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t m = get_mem_t_dim(a).y;
    size_t n = get_mem_t_dim(a).x;
    double* ad = current(a)(i,j);

    if((i+1)*m <= j*n) return;
    if(i*m >= (j+1)*n) return;
    for(size_t jj = j*n; jj < std::min((j+1)*n,(size_t)get_dim(a).x); jj++){
        if(i*m > jj) continue;
        if((i+1)*m <= jj) continue;
       *trace += ad[jj % m + (jj%n)*m];
    }
}

void transpose_c(pinned p_dense_matrix<double>& transposed, const p_dense_matrix<double>& original)
{
    int i = get_block_id(transposed).y;
    int j = get_block_id(transposed).x;
    double* td = current(transposed)(i,j);
    double* od = current(original)(j,i);

    for(size_t i = 0; i < get_mem_t_dim(original).y; ++i){
        for(size_t j=0; j < get_mem_t_dim(original).x; ++j){
            td[j+i*get_mem_t_dim(transposed).y] = od[i+j*get_mem_t_dim(original).y];
        }
    }
}

template <typename T>
void apply_writes_c(p_dense_matrix<T>& a)
{
    a.modifier = &a.modifiers.front();
    for(size_t k=0; k < a.modifier->size(); k++){
        size_t ii = ((*a.modifier)[k]).first.first;
        size_t jj = ((*a.modifier)[k]).first.second;
        size_t i = ii / get_mem_t_dim(a).y;
        size_t j = jj / get_mem_t_dim(a).x;
        if(!current(a).block(i,j)->available()) continue;
    
        T* ad = current(a)(i,j);
        size_t iii = ii % get_mem_t_dim(a).y;
        size_t jjj = jj % get_mem_t_dim(a).x;
        T value = *(T*)((*a.modifier)[k].second); 
        ad[iii + get_mem_t_dim(a).y * jjj] = value;
    }
    a.modifiers.pop();
    a.modifier = NULL;
}

void nullcut_c(pinned p_dense_matrix<double>& a, const size_t& num_rows, const size_t& num_cols)
{
    size_t i = get_block_id(a).y*get_mem_t_dim(a).y; 
    size_t j = get_block_id(a).x*get_mem_t_dim(a).x; 
    if((i+get_mem_t_dim(a).y <= num_rows) && (j+get_mem_t_dim(a).x <= num_cols)) return;

    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    for(size_t jj = 0; jj < get_mem_t_dim(a).x; jj++){
        if(j+jj < num_cols && (i+get_mem_t_dim(a).y > num_rows)) memset(&ad[jj*get_mem_t_dim(a).y+num_rows%get_mem_t_dim(a).y], 0, (get_mem_t_dim(a).y-num_rows%get_mem_t_dim(a).y)*sizeof(double));
        else if(j+jj >= num_cols) memset(&ad[jj*get_mem_t_dim(a).y], 0, get_mem_t_dim(a).y*sizeof(double));
    }
}
 
void validation_c(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, int*& bl) // C - See paper for Reference Dongara 
{ 
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x); 
    double* bd = current(b)(get_block_id(a).y, get_block_id(a).x); 
    double res(0.0); 
    double epsilon = std::numeric_limits<double>::epsilon();
    int position_x(0),position_y(0),position_xy(0); 

    for(int ii=0; ii < get_mem_t_dim(a).x; ++ii){ // C - the std::resize + cast dense to p_dense can add ghost element between lda + num_rows
        for(int jj=0; jj < get_mem_t_dim(a).y; ++jj){
            position_x = j*get_mem_t_dim(a).x+jj;
            position_y = i*get_mem_t_dim(a).x+ii;
            if(position_x < get_dim(a).x && position_y < get_dim(a).y){
                position_xy = ii*get_mem_t_dim(a).y+jj;
                res = (ad[position_xy]-bd[position_xy])/fabs(epsilon*bd[position_xy]); 
                if(res > 256){ // 16 is recommended by Dongara, 256 because the lapack give != runs after runs
                    printf("validation test failed in block %d %d, res %.16f matrix 1: %.16f matrix 2: %.16f \n", get_block_id(a).y, get_block_id(a).x, res, ad[i], bd[i]);
                *bl = 0; // test failed return 0 (bool false)
                }
            }
        }
    } 
}

// MKL LAPACK kernels
void svd_c(const p_dense_matrix<double>& a, int& m, int& n, p_dense_matrix<double>& u, p_dense_matrix<double>& vt, p_dense_matrix<double>& s)
{
/* Locals */
    int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
    int ldu = get_grid_dim(u).y*get_mem_t_dim(u).y;
    int ldvt = get_grid_dim(vt).y*get_mem_t_dim(vt).y;
    int info, lwork;
    double wkopt;
    double* work;
    assert(current(a).layout->get_list().size() != 0);
    current(a).solidify(current(a).layout->get_list());
    current(u).solidify(current(u).layout->get_list());
    current(vt).solidify(current(vt).layout->get_list());
    current(s).solidify(current(s).layout->get_list());
/* Query and allocate the optimal workspace */
    lwork = -1;
    dgesvd_( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(vt).data, &ldvt, &wkopt, &lwork, &info );
    lwork = (int)wkopt;
    work = (double*)malloc( lwork*sizeof(double) );
/* Compute SVD */
    dgesvd_( "S", "S", &m, &n, (double*)breakdown(a).data, &lda, (double*)breakdown(s).data, (double*)breakdown(u).data, &ldu, (double*)breakdown(vt).data, &ldvt, work, &lwork, &info );
/* Check for convergence */
    if( info > 0 ) {
        printf( "The algorithm computing SVD failed to converge.\n" );
        exit( 1 );
    }
    current(u).disperse(current(u).layout->get_list());
    current(vt).disperse(current(vt).layout->get_list());
    current(s).disperse(current(s).layout->get_list());
    free(work);
}

void syev_c(const p_dense_matrix<double>& a, int& m, p_dense_matrix<double>& w)
{
     int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
     int info, lwork = -1;

     double wkopt;
     double* work;
     assert(current(a).layout->get_list().size() != 0);
     current(a).solidify(current(a).layout->get_list());
     current(w).solidify(current(w).layout->get_list());
     
     dsyev_("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,&wkopt,&lwork,&info);
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     dsyev_("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,work,&lwork,&info);

     if( info > 0 ) {
         printf( "The algorithm computing SYEV failed to converge.\n" );
         exit( 1 );
     }

     current(a).disperse(current(a).layout->get_list());
     current(w).disperse(current(w).layout->get_list());
     free(work); 
}

void heev_c(const p_dense_matrix<double>& a, int& m, p_dense_matrix<double>& w)
{
     int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
     int info, lwork = -1;

     double wkopt;
     double* work;
     assert(current(a).layout->get_list().size() != 0);
     current(a).solidify(current(a).layout->get_list());
     current(w).solidify(current(w).layout->get_list());

     dsyev_("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,&wkopt,&lwork,&info);
     lwork = (int)wkopt;
     work = (double*)malloc( lwork*sizeof(double) );
     dsyev_("V","U",&m,(double*)breakdown(a).data,&lda,(double*)breakdown(w).data,work,&lwork,&info);

     if( info > 0 ) {
         printf( "The algorithm computing SYEV failed to converge.\n" );
         exit( 1 );
     }
    
     // First we reverse the eigenvalue, to be in agreement with the serial version ! 
     // The matrix is solidified, so we do not care on the workgroup representation
     double* dw = (double*)breakdown(w).data;
     double tempdbl;
     for (int i=0; i< static_cast<int>(m/2); i++){ 
         tempdbl = dw[i];
         dw[i] = dw[m-i-1];
         dw[m-i-1] = tempdbl;
     } 
     // Second we reverse the eigenvectors
     double* da = (double*)breakdown(a).data;
     double* tempcol = new double[lda]; 
     for (int i=0; i< static_cast<int>(m/2); ++i){ 
         memmove((void*)tempcol,(void*)&da[i*lda],lda*sizeof(double));
         memmove((void*)&da[i*lda],(void*)&da[(m-1-i)*lda],lda*sizeof(double));
         memmove((void*)&da[(m-1-i)*lda],(void*)tempcol,lda*sizeof(double));
     }
     delete[] tempcol; 
 
     current(a).disperse(current(a).layout->get_list());
     current(w).disperse(current(w).layout->get_list());
     free(work);
}

template<typename T>
void print_c(const p_dense_matrix<T>& a, int& m, int& n)
{
    int lda = get_grid_dim(a).y*get_mem_t_dim(a).y;
    assert(current(a).layout->get_list().size() != 0);

    current(a).solidify(current(a).layout->get_list());
    T* pa = (T*)breakdown(a).data;
    
    for(int i=0; i < m; ++i){
        for(int j=0; j < n; ++j){
            std::cout << *(pa+i+j*lda) << " ";                    
        }
        printf("\n");
    }

    current(a).disperse(current(a).layout->get_list());    
}
// C - For devs
void initv_c(pinned p_dense_matrix<double>&a, double const&v)
{
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    for(int i=0; i < get_mem_t_dim(a).x*get_mem_t_dim(a).y; ++i)
        ad[i] = v;
}

