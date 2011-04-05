#define NODE_COUNT 1
#include "mkl.h"

extern "C" {
    double  sqrt( double );
}

/*
 --- --- ---       --- --- ---       --- --- ---
| 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
| 0 | 1 | 2 |  x  | 0 | 1 | 2 |  =  | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
| 0 | 1 | 2 |     | 0 | 1 | 2 |     | 0 | 1 | 2 |
 --- --- ---       --- --- ---       --- --- ---
partial reduce?
*/

void gemm_c_kernel(pinned const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{
    int m   = get_group_t_dim(a).y;
    int n   = get_group_t_dim(b).x;
    int k   = get_group_t_dim(b).y;
    int lda = m;
    int ldb = k;
    int ldc = m;
    double alpha = 1.0; 
    double beta  = 0.0;
// a(i,j) => b(j,i) x a(z,j) where z : [0,m)
// current group of matrix a:
    int i = get_group_id(a).y;
    int j = get_group_id(a).x;
// taking (j,i) of b:
    double* bd = current(b)(j,i); // remote
// multiplying with column of a:
    for(int z = 0; z < get_grid_dim(a).y; z++){
        double* ad = current(a)(z,j);
        double* cd = reduced<'+'>(c)(z,i); // a(z,j) x b(j,i) => c(z,i)
        //printf("Performing DGEMM for %d,%d and %d,%d\n", z,j,j,i);
        dgemm("N","N", &m, &n, &k, &alpha, ad, &lda, bd, &ldb, &beta, cd, &ldc);
        //for(int ii=0; ii < m; ii++){
        //  for(int jj=0; jj < n; jj++)
        //  printf("%.2f	", cd[jj*ldc + ii]);
        //  printf("\n");
        //}
    }
}

void copy_c_kernel(p_dense_matrix<double>& ac, pinned const p_dense_matrix<double>& a)
{    
    int i = get_group_id(a).y;
    int j = get_group_id(a).x;
    double* a_elements  = current(a)(i,j);
    double* ac_elements = current(ac)(i,j);
    memcpy(ac_elements, a_elements, sizeof(double)*(get_group_dim(a)*get_item_dim(a)));
}

void remove_rows_c_kernel(pinned p_dense_matrix<double>& a, const size_t& i_mark, const size_t& k)
{
    typedef double T;
    
    double* ad    = NULL;
    double* ad_r  = NULL;
    double* ad_r0 = NULL;

    size_t i   = get_group_id(a).y;
    size_t j   = get_group_id(a).x;
    size_t lda = get_group_t_dim(a).y;

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
            for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                                           // memcpy from replacement block #1
                memcpy(&ad[lda*j + remains_u], &ad_r0[lda*j+lda-remains_l], sizeof(T)*remains_l);
        }else if(i >= group_i_mark){
            for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                                           // memcpy from replacement block #1
                memcpy(&ad[lda*j], &ad_r0[lda*j + (lda-remains)], sizeof(T)*remains);
        }
        for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                                               // memcpy from replacement block #2
            memcpy(&ad[lda*j + remains], &ad_r[lda*j], sizeof(T)*(lda-remains));
    }else{                                                                                             // get only one following block
        ambient::memoryfence();
        if(i == group_i_mark){
            if(remains_u + k < lda){
                for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                                       // first memmove inside block
                    memmove(&ad[lda*j + remains_u], &ad[lda*j + remains_u+k], sizeof(T)*(lda-remains_u-k));
                if(ad_r != NULL) for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                      // memcpy from replacement block
                    memcpy(&ad[lda*j + lda - k], &ad_r[lda*j], sizeof(T)*k);
            }else{
                if(ad_r != NULL) for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                      // memcpy from replacement block
                    memcpy(&ad[lda*j + remains_u], &ad_r[lda*j + lda-remains_l], sizeof(T)*(lda-remains_u));
            }
        }else if(i >= group_i_mark){
            for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                                           // first memmove inside block
                memmove(&ad[lda*j], &ad[lda*j + k_wo_blocks], sizeof(T)*(lda-k_wo_blocks) );
            if(ad_r != NULL) for(size_t j = 0; j < get_group_t_dim(a).x; ++j)                          // memcpy from replacement block
                memcpy(&ad[lda*j + lda-k_wo_blocks], &ad_r[lda*j], sizeof(T)*k_wo_blocks);
        }
    }
}

void remove_cols_c_kernel(pinned p_dense_matrix<double>& a, const size_t& j_mark, const size_t& k)
{
    typedef double T;

    double* ad    = NULL;
    double* ad_r  = NULL;
    double* ad_r0 = NULL;

    size_t i   = get_group_id(a).y;
    size_t j   = get_group_id(a).x;
    size_t lda = get_group_t_dim(a).y;
    size_t sda = get_group_t_dim(a).x;

    size_t remains_l = j_mark % sda;
    size_t remains_r = sda - (remains_l+k) % sda;
    size_t remains   = remains_l + remains_r;
    size_t shift     = __a_ceil(k / sda);
    size_t group_j_mark = j_mark / sda;
    size_t k_wo_blocks = std::min((2*sda-remains), k);

    if(j < group_j_mark) return;                                                                       // easy-out
    ad   = current(a)(i,j);
    if(j+shift < get_grid_dim(a).x) ad_r = current(a)(i,j+shift);

    if(remains < sda && (remains_l + k) > sda){                                                        // get two following blocks (j+shift-1;j+shift)
        if((j+shift-1) < get_grid_dim(a).x) ad_r0 = current(a)(i,j+shift-1);
        ambient::memoryfence();
        if(ad_r0 == NULL) return;                                                                      // out of matrix request
        if(j == group_j_mark){
            memcpy(&ad[lda*remains_l], &ad_r0[lda*(sda-remains_r)], sizeof(T)*lda*remains_r);          // memcpy from replacement block #1
        }else if(j >= group_j_mark){
            memcpy(ad, &ad_r0[lda*(sda-remains)], sizeof(T)*lda*remains);                              // memcpy from replacement block #1
        }
        memcpy(&ad[lda*remains], ad_r, sizeof(T)*lda*(sda-remains));                                   // memcpy from replacement block #2
    }else{                                                                                             // get only one following block
        ambient::memoryfence();
        if(j == group_j_mark){
            if(remains_l + k < sda){
                memmove(&ad[lda*remains_l], &ad[lda*(remains_l+k)], sizeof(T)*lda*(sda-remains_l-k));   // first memmove inside block
                if(ad_r != NULL) memcpy(&ad[lda*(sda-k)], ad_r, sizeof(T)*lda*k);                                       // memcpy from replacement block
            }else{
                if(ad_r != NULL) memcpy(&ad[lda*remains_l], &ad_r[lda*(sda-remains_r)], sizeof(T)*lda*(sda-remains_l)); // memcpy from replacement block
            }
        }else if(i >= group_j_mark){
            memmove(ad, &ad[lda*k_wo_blocks], sizeof(T)*lda*(sda-k_wo_blocks) );                        // first memmove inside block
            if(ad_r != NULL) memcpy(&ad[lda*(sda-k_wo_blocks)], ad_r, sizeof(T)*lda*k_wo_blocks);                       //  memcpy from replacement block
        }
    }
}

void resize_c_kernel(p_dense_matrix<double>& a, const size_t& rows, const size_t& cols)
{
    for(int i = 0; i < get_grid_dim(a).y; i++)
    if(current(a).group(i, get_grid_dim(a).x-1)->available()){
        size_t cutoff = get_grid_dim(a).x*get_group_t_dim(a).x - get_dim(a).x;
        if(cutoff > 0){
            double* ad = current(a)(i, get_grid_dim(a).x-1);
            memset(&ad[(get_group_t_dim(a).x - cutoff)*get_group_t_dim(a).y], 0, cutoff*sizeof(double));
        }
    }
    for(int j = 0; j < get_grid_dim(a).x; j++)
    if(current(a).group(get_grid_dim(a).y-1, j)->available()){
        double* ad = current(a)(get_grid_dim(a).y-1,j);
        size_t cutoff = get_grid_dim(a).y*get_group_t_dim(a).y - get_dim(a).y;
        size_t offset = get_group_t_dim(a).y - cutoff;
        if(cutoff > 0) for(int jj=0; jj < get_group_t_dim(a).x; jj++)
            memset(&ad[get_group_t_dim(a).y*jj+offset], 0, cutoff*sizeof(double));
    }
}

void sqrt_diagonal_c_kernel(pinned p_dense_matrix<double>& a){
    double* ad = current(a)(get_group_id(a).y, get_group_id(a).x);
    for(int i=0; i < get_group_t_dim(a).y; i++)
        ad[i] = sqrt(ad[i]);
}

void one_null_c_kernel(const p_dense_matrix<double>& a){ }

void two_null_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b){ }

void three_null_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, const p_dense_matrix<double>& c){ }

void add_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){
    double* ad = current(a)(get_group_id(c).y, get_group_id(c).x);
    double* cd = current(c)(get_group_id(c).y, get_group_id(c).x);
    printf("Executing add kernel for %d %d\n", get_group_id(c).y, get_group_id(c).x); 
    double* bd = current(b)(get_group_id(c).y, get_group_id(c).x);
    printf("Executed add kernel for %d %d\n", get_group_id(c).y, get_group_id(c).x); 
}
void sub_c_kernel(const p_dense_matrix<double>& a, const p_dense_matrix<double>& b, pinned p_dense_matrix<double>& c){}

void scale_c_kernel(const p_dense_matrix<double>& a, const double& alfa, pinned p_dense_matrix<double>& out){}
/**
void gemm_lhs_diagonal_c_kernel(const p_dense_matrix<double> & a_diag, pinned const p_dense_matrix<double>& b, p_dense_matrix<double>& c)
{

}*/

//void gemm_rhs_diagonal_c_kernel(pinned const p_dense_matrix<double> & a, const p_dense_matrix<double>& b_diag, p_dense_matrix<double>& c)
void gemm_rhs_diagonal_c_kernel(const p_dense_matrix<double> & a, const p_dense_matrix<double>& c)
{
	std::cout << " computation gemm diag " << std::endl;  
 //  int j = get_group_id(a).x*get_group_t_dim(a).x;
/*
    int size = get_group_t_dim(a).y;
    int ONE = 1;
    double* ad = current(a)(get_group_id(a).y, get_group_id(a).x);
    double* cd = current(c)(get_group_id(a).y, get_group_id(a).x);

    memset(cd, 0, get_group_t_dim(c).x*get_group_t_dim(c).y*sizeof(double));
    for(int jj = 0 ; jj < get_group_t_dim(a).x ; jj++)
    {
	 double * alpha = current(b_diag)(j/get_group_t_dim(b_diag).y);
	 daxpy(&size,&alpha[j%get_group_t_dim(b_diag).y],&ad[jj*get_group_t_dim(a).y], &ONE,&cd[jj*get_group_t_dim(c).y],&ONE);
    } 
*/
}



void single_integer_c_kernel(int*& input){
    zout << "single integer kernel: output is " << input[2] << "\n";
}
