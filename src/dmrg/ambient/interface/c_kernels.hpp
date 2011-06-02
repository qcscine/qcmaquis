#define NODE_COUNT 1

//#include "mkl.h"

extern "C" {
    double sqrt(double);
    double fabs(double);
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
    //printf("copy_c\n");
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    if(i >= get_grid_dim(ac).y || j >= get_grid_dim(ac).x) return;
    double* a_elements  = current(a)(i,j);
    double* ac_elements = current(ac)(i,j);
    memcpy(ac_elements, a_elements, sizeof(double)*get_mem_t_dim(a).y*get_mem_t_dim(a).x);
}

/** Merge Sort algo from John von Neumann, divide and conquer.
http://en.literateprograms.org/Merge_sort_(C_Plus_Plus) **/
template<class T>
void merge_sort(T begin, T end)
{
    size_t size(end-begin);
    if(size < 2) return; // end recurrence 
   
    T begin_right = begin+size/2; 
    merge_sort(begin, begin_right);
    merge_sort(begin_right, end);
    merge(begin, begin_right, end);
}

void associated_sort_c(pinned p_dense_matrix<double>& a)
{    
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t mem_t_dim_y = get_mem_t_dim(a).y;
    size_t grid_dim_y = get_grid_dim(a).y;
    size_t offset = 0;

    // pointers for the first merge sort
    double* ad_begin = current(a)(i,j);
    double* ad_end   = ad_begin + mem_t_dim_y;
    // pointers for the odd/even sort
    double* ad;
    double* ad_oe;
    /** von Neumann sort on every workblock **/
    merge_sort(ad_begin, ad_end);
    /** Final sort between work group, A new Parallel Sorting Algo based on odd-even Mergesort,
    15th euromicro Int. Conf. on Parallel, Distributed-Based Processing **/    

    /*for(size_t ii = 0; ii < (grid_dim_y/2);ii++){
        if((i%2) == 0){ // even
            if(i<(grid_dim_y-1)){ // to avoid to request a none exsiting block
                 ad = current(a)(i,j);
                 //ad_oe = modified(a)(i,j);
                 ad_oe = current(a)(i+1,j);
                 //ad_oe = modified(a)(i+1,j);
                 mpi_merge(ad,ad_oe,mem_t_dim_y);
             } 
        }else{ // odd
             if(i<(grid_dim_y-1)){ // to avoid to request a none exsiting block
                 ad = current(a)(i,j);
                 //ad_oe = modified(a)(i,j);
                 ad_oe = current(a)(i+1,j);
                 //ad_oe = modified(a)(i+1,j);
                 mpi_merge(ad,ad_oe,mem_t_dim_y);
             }
        }
    }*/
}
//////////////////////////////////////////////////////
template<class T, class VT>
void insert(T begin, T end, VT &v)
{
    using ::std::swap;
    while(begin+1!=end && *(begin+1)<v){
        swap(*begin, *(begin+1));
        begin++;
    }
    swap(*begin,v);
}

template<class T>
void merge(T begin, T begin_right, T end)
{
    for(;begin<begin_right;begin++){// need contiguous mem
        if(*begin_right < *begin){
            typename std::iterator_traits<T>::value_type v(0.0);
            using ::std::swap;
            swap(v, *begin);
            swap(*begin, *begin_right);
            insert(begin_right, end, v);
        }
    }
}

/** we need the mpi_merge when the data is not contiguous **/
template<class T>
void mpi_merge(T left, T right, size_t size)
{
   for(size_t jj=0; jj<size; jj++){// we can not used directly merge, because, we have not contiguous mem  
       if(*right < *(left+jj)){
           typename std::iterator_traits<T>::value_type v(0.0);
           using ::std::swap;
           swap(v,*(left+jj));
           swap(*(left+jj),*(right));
           insert(right,(right+size),v); // pointer sum
       }
   }
}

void associated_sort_e_c(pinned p_dense_matrix<double>& a)
{
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t mem_t_dim_y = get_mem_t_dim(a).y;
    size_t grid_dim_y = get_grid_dim(a).y;

    if((i%2) == 0){
        if(i<(grid_dim_y-1)){ // to avoid to request a none exsiting block
             double* ad = current(a)(i,j);
             double* ad_oe = current(a)(i+1,j);
             mpi_merge(ad,ad_oe,mem_t_dim_y);
        }
    }
}

void associated_sort_o_c(pinned p_dense_matrix<double>& a)
{
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t mem_t_dim_y = get_mem_t_dim(a).y;
    size_t grid_dim_y = get_grid_dim(a).y;

    if((i%2) != 0){
        if(i<(grid_dim_y-1)){ // to avoid to request a none exsiting block
             double* ad = current(a)(i,j);
             double* ad_oe = current(a)(i+1,j);
             mpi_merge(ad,ad_oe,mem_t_dim_y);
        }
    }
}

void associated_reverse_c(p_dense_matrix<double>& a, const size_t& num_rows)
{   
    double rs = 0; // intermediate for reverse inside one block of mem
    size_t iir = 0;
 
    size_t grid_dim_y = get_grid_dim(a).y;

    size_t size = get_mem_t_dim(a).y;
    size_t size_total = get_mem_t_dim(a).y*grid_dim_y; 
    //first we reverse all elements of block dim
    for(size_t jj = 0 ; jj < grid_dim_y  ; jj++ ) {
        double* ad = current(a)(jj,0);

        #if defined (INTEL_COMPILER)
        #pragma ivdep // should be ok
        #endif
        for(size_t ii = 0; ii < size/2; ii++){  
            iir = size - ii - 1; // for SIMD/Stream
            rs  = ad[ii];
            ad[ii]  = ad[iir];
            ad[iir] = rs;
        }  
    }
    
    double* pa = (double*)malloc(size*sizeof(double));
    int jjr;
    // second we reverse the blocks of mem 
    for(size_t jj = 0; jj < grid_dim_y/2; jj++){
        jjr = grid_dim_y - jj - 1;
        memcpy((void*)pa, (void*)current(a)(jj,0), size*sizeof(double));   
        memcpy((void*)current(a)(jj,0), (void*)(current(a)(jjr,0)), size*sizeof(double));
        memcpy((void*)current(a)(jjr,0), (void*)pa, size*sizeof(double));
    } 

    // finalization due to a possible offset
    size_t markmem = size_total - num_rows;
    size_t offset  = size - markmem; 
    if(markmem != 0){ // we have an offset   
        for(size_t jj = 0; jj < grid_dim_y; jj++){ // first loop, offset intra block
           double* ad = current(a)(jj,0);
            memmove((void*)ad, (void*)&ad[markmem], offset*sizeof(double)); // intrablock shift
            if(jj != grid_dim_y - 1){
                memmove((void*)&ad[offset], (void*)current(a)(jj+1,0), markmem*sizeof(double)); // extrablock shift
            }else{
                memset((void*)&ad[markmem], 0, offset*sizeof(double)); // fill up with 0 until the end
            }
        } 
    }
    
    free((void*)pa);
}

void move_offset_c(pinned p_dense_matrix<double>& a)
{
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t offset = get_grid_dim(a).y*get_mem_t_dim(a).y - get_dim(a).y; 
    if(offset == 0) return; // no offset 
    size_t mem_y = get_mem_t_dim(a).y - offset;
    size_t mem_t_dim_y = get_mem_t_dim(a).y;
    double* adp;
    double* ad = current(a)(i,j);

    if(i < (get_grid_dim(a).y-1)){
        adp = current(a)((i+1),j);
        memmove((void*)ad,(void*)(ad+offset),mem_y*sizeof(double));
        memcpy((void*)(ad+mem_y),(void*)(adp),offset*sizeof(double));
    }
    if(i == get_grid_dim(a).y-1){ //final offset
        memmove((void*)ad,(void*)(ad+offset),mem_y*sizeof(double));
    }
}

void associated_find_if_c(pinned p_dense_matrix<double>& a, const double& value, size_t*& out_value)
{ // only single process is supported (outmost assign)
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    double* ad = current(a)(i,j);
   
    size_t ii=0;
    while(ad[ii]>value && ii < get_mem_t_dim(a).y) {
        *out_value = ii+i*get_mem_t_dim(a).y;
        ii++;
    }
}
 
void associated_accumulate_c(pinned p_dense_matrix<double>& a, const size_t*& begin, double*& out_value)
{
    /** to check **/
    size_t i = get_block_id(a).y;
    size_t j = get_block_id(a).x;
    size_t mem_t_dim_y = get_mem_t_dim(a).y;
    double* ad = current(a)(i,j);
    double sum(0.0);    

    for(size_t ii=0; ii<mem_t_dim_y;ii++){
        if(i*mem_t_dim_y+ii > *begin){
            sum += ad[ii];
        } 
    }
    *out_value = sum;
}

void associated_max_c(pinned p_dense_matrix<double>&a, const double& evalscut, const size_t& Mmax, double* & out_value)
{
    size_t i = Mmax/get_mem_t_dim(a).y; // we are looking where is the corresponding element to S[Mmax]
    size_t j = 0;
    size_t ii = Mmax%get_mem_t_dim(a).y;   
    double* ad = current(a)(i,j);

    if(ad[ii]>evalscut){
        *out_value=ad[ii];
    }else{
        *out_value=evalscut;
    }
}

void variable_free_c(void*& a){ free(a); }

void remove_rows_c(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    typedef double T;
    
    double* ad    = NULL;
    double* ad_r  = NULL;
    double* ad_r0 = NULL;

    size_t i   = get_block_id(a).y;
    size_t j   = get_block_id(a).x;
    size_t lda = get_mem_t_dim(a).y;

    size_t remains_u = i_mark % lda;
    size_t remains_l = lda - (remains_u+k) % lda;
    size_t remains   = remains_u + remains_l;
    size_t shift     = __a_ceil(k / lda);
    size_t group_i_mark = i_mark / lda;
    size_t k_wo_blocks = std::min((2*lda-remains), k);
    if(i < group_i_mark) return;                                                                       // easy-out
    ad   = current(a)(i,j);
    if(i+shift < get_grid_dim(a).y) ad_r = current(a)(i+shift,j);
 
    if(remains < lda && (remains_u + k) > lda){                                                        // get two following blocks (i+shift-1;i+shift)
        if((i+shift-1) < get_grid_dim(a).y) ad_r0 = current(a)(i+shift-1,j);
        ambient::memoryfence();
        if(ad_r0 == NULL) return;                                                                      // out of matrix request
        if(i == group_i_mark){
            for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                                             // memcpy from replacement block #1
                memcpy(&ad[lda*j + remains_u], &ad_r0[lda*j+lda-remains_l], sizeof(T)*remains_l);
        }else if(i >= group_i_mark){
            for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                                             // memcpy from replacement block #1
                memcpy(&ad[lda*j], &ad_r0[lda*j + (lda-remains)], sizeof(T)*remains);
        }
        for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                                                 // memcpy from replacement block #2
            memcpy(&ad[lda*j + remains], &ad_r[lda*j], sizeof(T)*(lda-remains));
    }else{                                                                                             // get only one following block
        ambient::memoryfence();
        if(i == group_i_mark){
            if(remains_u + k < lda){
                for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                                         // first memmove inside block
                    memmove(&ad[lda*j + remains_u], &ad[lda*j + remains_u+k], sizeof(T)*(lda-remains_u-k));
                if(ad_r != NULL) for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                        // memcpy from replacement block
                    memcpy(&ad[lda*j + lda - k], &ad_r[lda*j], sizeof(T)*k);
            }else{
                if(ad_r != NULL) for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                        // memcpy from replacement block
                    memcpy(&ad[lda*j + remains_u], &ad_r[lda*j + lda-remains_l], sizeof(T)*(lda-remains_u));
            }
        }else if(i >= group_i_mark){
            for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                                             // first memmove inside block
                memmove(&ad[lda*j], &ad[lda*j + k_wo_blocks], sizeof(T)*(lda-k_wo_blocks));
            if(ad_r != NULL) for(size_t j = 0; j < get_mem_t_dim(a).x; ++j)                            // memcpy from replacement block
                memcpy(&ad[lda*j + lda-k_wo_blocks], &ad_r[lda*j], sizeof(T)*k_wo_blocks);
        }
    }
}

void remove_cols_c(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{
    typedef double T;

    double* ad    = NULL;
    double* ad_r  = NULL;
    double* ad_r0 = NULL;

    size_t i   = get_block_id(a).y;
    size_t j   = get_block_id(a).x;
    size_t lda = get_mem_t_dim(a).y;
    size_t sda = get_mem_t_dim(a).x;

    size_t remains_l = j_mark % sda;
    size_t remains_r = sda - (remains_l+k) % sda;
    size_t remains   = remains_l + remains_r;
    size_t shift     = __a_ceil(k / sda);
    size_t group_j_mark = j_mark / sda;
    size_t k_wo_blocks = std::min((2*sda-remains), k);

    if(j < group_j_mark) return;                                                                                        // easy-out
    ad   = current(a)(i,j);                                                                                        
    if(j+shift < get_grid_dim(a).x) ad_r = current(a)(i,j+shift);                                                  
                                                                                                                   
    if(remains < sda && (remains_l + k) > sda){                                                                         // get two following blocks (j+shift-1;j+shift)
        if((j+shift-1) < get_grid_dim(a).x) ad_r0 = current(a)(i,j+shift-1);                                       
        ambient::memoryfence();                                                                                    
        if(ad_r0 == NULL) return;                                                                                       // out of matrix request
        if(j == group_j_mark){                                                                                     
            memcpy(&ad[lda*remains_l], &ad_r0[lda*(sda-remains_r)], sizeof(T)*lda*remains_r);                           // memcpy from replacement block #1
        }else if(j >= group_j_mark){                                                                               
            memcpy(ad, &ad_r0[lda*(sda-remains)], sizeof(T)*lda*remains);                                               // memcpy from replacement block #1
        }                                                                                                          
        memcpy(&ad[lda*remains], ad_r, sizeof(T)*lda*(sda-remains));                                                    // memcpy from replacement block #2
    }else{                                                                                                              // get only one following block
        ambient::memoryfence();
        if(j == group_j_mark){
            if(remains_l + k < sda){
                memmove(&ad[lda*remains_l], &ad[lda*(remains_l+k)], sizeof(T)*lda*(sda-remains_l-k));                   // first memmove inside block
                if(ad_r != NULL) memcpy(&ad[lda*(sda-k)], ad_r, sizeof(T)*lda*k);                                       // memcpy from replacement block
            }else{
                if(ad_r != NULL) memcpy(&ad[lda*remains_l], &ad_r[lda*(sda-remains_r)], sizeof(T)*lda*(sda-remains_l)); // memcpy from replacement block
            }
        }else if(i >= group_j_mark){
            memmove(ad, &ad[lda*k_wo_blocks], sizeof(T)*lda*(sda-k_wo_blocks) );                                        // first memmove inside block
            if(ad_r != NULL) memcpy(&ad[lda*(sda-k_wo_blocks)], ad_r, sizeof(T)*lda*k_wo_blocks);                       //  memcpy from replacement block
        }
    }
}

void touch_c(p_dense_matrix<double>& a){ }

void resize_c(p_dense_matrix<double>& a, const size_t& rows, const size_t& cols)
{
    for(int i = 0; i < get_grid_dim(a).y; i++)
    if(current(a).block(i, get_grid_dim(a).x-1)->available()){
        size_t cutoff = get_grid_dim(a).x*get_mem_t_dim(a).x - get_dim(a).x;
        if(cutoff > 0){
            double* ad = current(a)(i, get_grid_dim(a).x-1);
            memset(&ad[(get_mem_t_dim(a).x - cutoff)*get_mem_t_dim(a).y], 0, cutoff*sizeof(double));
        }
    }
    for(int j = 0; j < get_grid_dim(a).x; j++)
    if(current(a).block(get_grid_dim(a).y-1, j)->available()){
        double* ad = current(a)(get_grid_dim(a).y-1,j);
        size_t cutoff = get_grid_dim(a).y*get_mem_t_dim(a).y - get_dim(a).y;
        size_t offset = get_mem_t_dim(a).y - cutoff;
        if(cutoff > 0) for(int jj=0; jj < get_mem_t_dim(a).x; jj++)
            memset(&ad[get_mem_t_dim(a).y*jj+offset], 0, cutoff*sizeof(double));
    }
}

void sqrt_diagonal_c(pinned p_dense_matrix<double>& a)
{
    double* ad = current(a)(get_block_id(a).y, get_block_id(a).x);
    for(int i=0; i < get_mem_t_dim(a).y; i++)
    ad[i] = sqrt(ad[i]);
}

template<typename T>
void __a_memcpy(T& dest, dim2 dest_p, const T& src, dim2 src_p, dim2 size)
{
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

    for(size_t j = startj; j < limj; j++){
        size_t sj = dj + j - dest_p.x + src_p.x;
        size_t sii = si % get_mem_t_dim(src).y;
        size_t sjj = sj % get_mem_t_dim(src).x;
        size_t w = limi - starti;
        size_t i = starti;
        for(int k = 0; k < src_blocks_i; k++){
            typename T::value_type* sd = current(src)(si / get_mem_t_dim(src).y + k, sj / get_mem_t_dim(src).x);
            memcpy(&dd[j*get_mem_t_dim(dest).y + i],
                   &sd[sjj*get_mem_t_dim(src).y+sii],
                   std::min(get_mem_t_dim(src).y-sii, w)*sizeof(typename T::value_type));
            w -= get_mem_t_dim(src).y-sii;
            i += get_mem_t_dim(src).y-sii;
            sii = 0;
        }
    }
}

void copy_after_c(pinned p_dense_matrix<double>& ac, const size_t& pos, const p_dense_matrix<double>& a)
{
    __a_memcpy(ac, dim2(0,pos), a, dim2(0,0), dim2(1,get_dim(a).y));
}

void reshape_l2r_c(const p_dense_matrix<double>& left, pinned p_dense_matrix<double>& right,
                          const size_t& left_offset, const size_t& right_offset, 
                          const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    //printf("reshape_l2r_c\n");
    for(size_t ss = 0; ss < sdim; ++ss)
        __a_memcpy(right, dim2(ss*rdim + right_offset,0), 
                   left,  dim2(0, ss*ldim + left_offset), 
                   dim2( rdim, ldim ));
}

void reshape_r2l_c(pinned p_dense_matrix<double>& left, const p_dense_matrix<double>& right,
                          const size_t& left_offset, const size_t& right_offset, 
                          const size_t& sdim, const size_t& ldim, const size_t& rdim)
{
    for(size_t ss = 0; ss < sdim; ++ss)
        __a_memcpy(left,  dim2(0, ss*ldim + left_offset), 
                   right, dim2(ss*rdim + right_offset,0), 
                   dim2( rdim, ldim ));
}

template<typename T>
void __a_add_scaled(T& dest, dim2 dest_p, const T& src, dim2 src_p, typename T::value_type alfa, dim2 size)
{
    size_t starti, startj, limi, limj;
    size_t di = get_block_id(dest).y * get_mem_t_dim(dest).y;
    size_t dj = get_block_id(dest).x * get_mem_t_dim(dest).x;

    assert(get_grid_dim(dest).x*get_mem_t_dim(dest).x - dest_p.x >= size.x);
    if(get_grid_dim(dest).y*get_mem_t_dim(dest).y - dest_p.y < size.y) printf("%d .. %d; %d vs %d\n", (int)(get_grid_dim(dest).y*get_mem_t_dim(dest).y - dest_p.y), (int)size.y, (int)dest.num_rows(), (int)get_grid_dim(dest).y*get_mem_t_dim(dest).y );
    assert(get_grid_dim(dest).y*get_mem_t_dim(dest).y - dest_p.y >= size.y);
    assert(get_grid_dim(src).x*get_mem_t_dim(src).x - src_p.x >= size.x);
    if(get_grid_dim(src).y*get_mem_t_dim(src).y - src_p.y < size.y) printf("%d .. %d\n", (int)(get_grid_dim(src).y*get_mem_t_dim(src).y - src_p.y), (int)size.y );
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
    size_t src_blocks_i = 0;
    int num_src_blocks = limi-starti-get_mem_t_dim(src).y+sii;
    if(num_src_blocks > 0) src_blocks_i = __a_ceil( num_src_blocks / get_mem_t_dim(src).y ) + 1;
// let's exhaust first src block
    typename T::value_type* dd = current(dest)(get_block_id(dest).y, get_block_id(dest).x);

    for(size_t j = startj; j < limj; j++){
        size_t sj = dj + j - dest_p.x + src_p.x;
        size_t sii = si % get_mem_t_dim(src).y;
        size_t sjj = sj % get_mem_t_dim(src).x;
        size_t w = limi - starti;
        size_t i = starti;
        for(int k = 0; k < src_blocks_i; k++){
            typename T::value_type* sd = current(src)(si / get_mem_t_dim(src).y + k, sj / get_mem_t_dim(src).x);
            for(int z = 0; z < std::min(get_mem_t_dim(src).y-sii, w); z++)
                dd[j*get_mem_t_dim(dest).y+i + z] += sd[sjj*get_mem_t_dim(src).y+sii + z]*alfa;
            w -= get_mem_t_dim(src).y-sii;
            i += get_mem_t_dim(src).y-sii;
            sii = 0;
        }
    }
}

template <typename T>
void rb_tensor_mpo_c(pinned p_dense_matrix<T>& out, const p_dense_matrix<T>& in, const p_dense_matrix<T>& alfa,
                          const size_t& out_offset, const size_t& in_offset, 
                          const size_t& sdim1, const size_t& sdim2, const size_t& ldim, const size_t& rdim)
{
    //printf("rb_tensor_mpo\n");
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
    //printf("rb_tensor_mpo\n");
    for(size_t ss1 = 0; ss1 < sdim1; ++ss1)
        for(size_t ss2 = 0; ss2 < sdim2; ++ss2){
            T* alfad = current(alfa)(ss1/get_mem_t_dim(alfa).y, ss2/get_mem_t_dim(alfa).x);
            T  alfa_t = alfad[ss1%get_mem_t_dim(alfa).y + get_mem_t_dim(alfa).y*(ss2%get_mem_t_dim(alfa).x)];
            __a_add_scaled(out, dim2(0, out_offset + ss2*ldim),
                           in,  dim2(0, in_offset + ss1*ldim),
                           alfa_t, dim2(rdim, ldim));
        }
}

void scalar_norm_c(pinned const p_dense_matrix<double>& a, p_dense_matrix<double>& norm)
{
    double summ = 0;
    int i = get_block_id(a).y;
    int j = get_block_id(a).x;
    double* ad = current(a)(i,j);
    for(size_t ii=0; ii < get_mem_t_dim(a).x*get_mem_t_dim(a).y; ii++)
        summ += ad[ii]*ad[ii];

    double* ret = reduced<'+'>(norm)(0,0);
    *ret = summ;
}

void atomic_add_c(p_dense_matrix<double>& a, const p_dense_matrix<double>& b)
{
    double* ad = current(a)(0,0);
    double* bd = current(b)(0,0);
    *ad += *bd;
}

void add_c(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c)
{
    double* ad = current(a)(get_block_id(c).y, get_block_id(c).x);
    double* bd = current(b)(get_block_id(c).y, get_block_id(c).x);
    double* cd = current(c)(get_block_id(c).y, get_block_id(c).x);
    for(int i=0; i < get_mem_t_dim(c).x*get_mem_t_dim(c).y; i++)
    cd[i] = ad[i] + bd[i];
}

void sub_c(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c)
{
    double* ad = current(a)(get_block_id(c).y, get_block_id(c).x);
    double* bd = current(b)(get_block_id(c).y, get_block_id(c).x);
    double* cd = current(c)(get_block_id(c).y, get_block_id(c).x);
    for(int i=0; i < get_mem_t_dim(c).x*get_mem_t_dim(c).y; i++)
    cd[i] = ad[i] - bd[i];
}

void scale_c(const p_dense_matrix<double>& m, const double& t, pinned p_dense_matrix<double>& out)
{
    double* md   = current(m)(get_block_id(out).y, get_block_id(out).x);
    double* outd = current(out)(get_block_id(out).y, get_block_id(out).x);
    for(int i=0; i < get_mem_t_dim(out).x*get_mem_t_dim(out).y; i++)
    outd[i] = md[i]*t;
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
 
void one_null_c(const p_dense_matrix<double>& a){ }
void two_null_c(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b){ }
void three_null_c(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, const p_dense_matrix<double>& c){ }

/** 
Validation kernel source:  
 @article{Dongarra:1990:SLB:77626.79170, 
 author = {Dongarra, J. J. and Du Croz, Jeremy and Hammarling, Sven and Duff, I. S.}, 
 title = {A set of level 3 basic linear algebra subprograms}, 
 journal = {ACM Trans. Math. Softw.}, 
 volume = {16}, 
 issue = {1}, 
 month = {March}, 
 year = {1990}, 
 issn = {0098-3500}, 
 pages = {1--17}, 
 numpages = {17}, 
 url = {http://doi.acm.org/10.1145/77626.79170}, 
 doi = {http://doi.acm.org/10.1145/77626.79170}, 
 acmid = {79170}, 
 publisher = {ACM}, 
 address = {New York, NY, USA}, 
}  
*/ 
void validation_c(pinned const p_dense_matrix<double>& a_ambient, const p_dense_matrix<double>& b_scalapack) 
{ 
    double* ad = current(a_ambient)(get_block_id(a_ambient).y, get_block_id(a_ambient).x); 
    double* bd = current(b_scalapack)(get_block_id(a_ambient).y, get_block_id(a_ambient).x); 

    double res = 0; 
    double epsilon = 0.0000000000000001; 
    for(int i=0; i < get_mem_t_dim(a_ambient).x*get_mem_t_dim(a_ambient).y; i++) 
    { 
        res = (fabs(ad[i]-bd[i]))/fabs(epsilon*bd[i]); 
        if(res > 16){ // 16 is recommended by Dongara,  
             printf("validation test failed in block %d %d, res %.10f Ambient: %.10f Scala: %.10f \n", get_block_id(a_ambient).y, get_block_id(a_ambient).x, res, ad[i], bd[i]);
             //exit(-1);
        }                
    } 
}

void associated_validation_c(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b) 
{// two validations kernels is stupid, untill the implementation of the << >> operator inside the l_kernel 
    int j = get_block_id(a).x;
    int i = get_block_id(a).y;

    double* ad = current(a)(i,j);
    double* bd = current(b)(i,j);

    double res = 0; 
    double epsilon = 0.0000000000000001; 
  
    for(int jj=0; jj < get_mem_t_dim(a).y; jj++) 
    { 
        res = (fabs(ad[jj]-bd[jj]))/fabs(epsilon*bd[jj]); 
        if(res > 16){ // 16 is recommended by dongara,  
             printf("validation test failed in block %d %d, res %.10f Ambient: %.10f Scala: %.10f \n", get_block_id(a).y, get_block_id(a).x, res, ad[jj], bd[jj]);
             exit(-1);
        }                
    } 
} 
