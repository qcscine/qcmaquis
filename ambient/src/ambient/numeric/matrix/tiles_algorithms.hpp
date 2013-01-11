#ifndef __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__
#define __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__

#include "ambient/numeric/matrix/tiles.h"

#define value_type      typename tiles<Matrix>::value_type
#define size_type       typename tiles<Matrix>::size_type
#define real_type       typename tiles<Matrix>::real_type
#define scalar_type     typename tiles<Matrix>::scalar_type
#define difference_type typename tiles<Matrix>::difference_type

inline size_t __a_mod_classic(size_t size, size_t tile){
    size_t mask[2] = {(size & (tile-1)), tile}; 
    return mask[!mask[0]];
}

inline size_t __a_mod(size_t size, size_t tile){
    size_t m = size % tile;
    if(m == 0) m = tile;
    return m;
}

template<typename T>
inline void __a_reduce(std::vector<T*>& seq){
    if(seq.size() == 1) return;
    for(int stride = 1; stride < seq.size(); stride *= 2)
        for(int k = stride; k < seq.size(); k += stride*2){
            *seq[k-stride] += *seq[k];
        }
}

class cross_iterator {
public:
    cross_iterator(size_t first, size_t second, size_t size) 
    : first(first), second(second), lim(first+size){
        measure_step();
    }
    void operator++ (){
        first  += step;
        second += step;
        measure_step();
    }
    bool end(){
        return (first >= lim);
    }
    void measure_step(){
        step = std::min(std::min((AMBIENT_IB*__a_ceil((first+1)/AMBIENT_IB) - first), 
                                 (AMBIENT_IB*__a_ceil((second+1)/AMBIENT_IB) - second)),
                                 (lim-first));
    }
    size_t first;
    size_t second;
    size_t step;
    size_t lim;
};

namespace ambient { namespace numeric {

    template<class Matrix>
    bool is_hermitian(const tiles<Matrix>& a){
        return false;
    }

    template<class Matrix>
    inline void split_d(const tiles<Matrix>& a){
        if(a.data.size() != 1) return;
        if(a.mt == 1 && a.nt == 1) return;

        std::vector<Matrix*> split;
        int tailm = __a_mod(a.rows, AMBIENT_IB);
        int tailn = __a_mod(a.cols, AMBIENT_IB);
        split.reserve(a.mt*a.nt);
        for(int j = 1; j < a.nt; j++){
            for(int i = 1; i < a.mt; i++) 
                split.push_back(new Matrix(AMBIENT_IB, AMBIENT_IB));
            split.push_back(new Matrix(tailm, AMBIENT_IB));
        }
        for(int i = 1; i < a.mt; i++) 
            split.push_back(new Matrix(AMBIENT_IB, tailn));
        split.push_back(new Matrix(tailm, tailn));

        const Matrix* src = a.data[0];

        if(!a[0].impl->weak())
        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                Matrix& dst = *split[i+j*a.mt];
                copy_block(*src, i*AMBIENT_IB, j*AMBIENT_IB, dst, 0, 0, dst.num_rows(), dst.num_cols());
            }
        }

        delete src;
        std::swap(const_cast<std::vector<Matrix*>&>(a.data), split);
    }

    template<class Matrix>
    inline void split_d(const tiles<transpose_view<Matrix> >& a){
        if(a.data.size() == 1 && (__a_ceil(a.rows/AMBIENT_IB) != 1 || __a_ceil(a.cols/AMBIENT_IB) != 1)) printf("We have to split but we don't :(\n");
    }

    template<typename T>
    inline void split_d(const tiles<diagonal_matrix<T> >& a){
        if(a.data.size() != 1) return;
        if(a.nt == 1) return;

        std::vector<diagonal_matrix<T> *> split;
        int tailm = __a_mod(a.size, AMBIENT_IB);
        split.reserve(a.nt);
        for(int i = 1; i < a.nt; i++) 
            split.push_back(new diagonal_matrix<T>(AMBIENT_IB));
        split.push_back(new diagonal_matrix<T>(tailm));

        const diagonal_matrix<T>* src = a.data[0];

        if(!a[0].get_data().impl->weak())
        for(int i = 0; i < a.nt; i++){
            diagonal_matrix<T>& dst = *split[i];
            copy_block(src->get_data(), i*AMBIENT_IB, 0, dst.get_data(), 0, 0, dst.num_rows(), 1);
        }

        delete src;
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), split);
    }
    
    template<class Matrix>
    inline Matrix* merge(const tiles<subset_view<Matrix> >& a){
        Matrix* m = new Matrix(a.rows, a.cols);

        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                const Matrix* src = &a.tile(i,j);
                copy_block(*src, 0, 0, *m, i*AMBIENT_IB, j*AMBIENT_IB, src->num_rows(), src->num_cols());
            }
        }
        return m;
    }

    // {{{ normal merge / split
    
    template<class Matrix>
    inline void merge(const tiles<Matrix>& a){
        if(a.data.size() == 1) return;

        std::vector<Matrix*> m; 
        m.push_back(new Matrix(a.rows, a.cols));

        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                const Matrix* src = &a.tile(i,j);
                if(!src->impl->weak())
                copy_block(*src, 0, 0, *m[0], i*AMBIENT_IB, j*AMBIENT_IB, src->num_rows(), src->num_cols());
                delete src;
            }
        }
        std::swap(const_cast<std::vector<Matrix*>&>(a.data), m);
    }

    template<class Matrix>
    inline void split(const tiles<Matrix>& a){
        return;
        if(a.data.size() != 1) return;
        if(a.mt == 1 && a.nt == 1) return;

        std::vector<Matrix*> split;
        int tailm = __a_mod(a.rows, AMBIENT_IB);
        int tailn = __a_mod(a.cols, AMBIENT_IB);
        split.reserve(a.mt*a.nt);
        for(int j = 1; j < a.nt; j++){
            for(int i = 1; i < a.mt; i++) 
                split.push_back(new Matrix(AMBIENT_IB, AMBIENT_IB));
            split.push_back(new Matrix(tailm, AMBIENT_IB));
        }
        for(int i = 1; i < a.mt; i++) 
            split.push_back(new Matrix(AMBIENT_IB, tailn));
        split.push_back(new Matrix(tailm, tailn));

        const Matrix* src = a.data[0];

        if(!a[0].impl->weak())
        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                Matrix& dst = *split[i+j*a.mt];
                copy_block(*src, i*AMBIENT_IB, j*AMBIENT_IB, dst, 0, 0, dst.num_rows(), dst.num_cols());
            }
        }

        delete src;
        std::swap(const_cast<std::vector<Matrix*>&>(a.data), split);
    }
    // }}}
    // {{{ diagonal merge/split
    template<typename T>
    inline void merge(const tiles<diagonal_matrix<T> >& a){
        if(a.data.size() == 1) return;

        std::vector<diagonal_matrix<T> *> m; 
        m.push_back(new diagonal_matrix<T>(a.size));

        for(int i = 0; i < a.nt; i++){
            const diagonal_matrix<T>* src = a.data[i];
            if(!src->get_data().impl->weak())
            copy_block(src->get_data(), 0, 0, m[0]->get_data(), i*AMBIENT_IB, 0, src->num_rows(), 1);
            delete src;
        }
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), m);
    }

    template<typename T>
    inline void split(const tiles<diagonal_matrix<T> >& a){
        return;
        if(a.data.size() != 1) return;
        if(a.nt == 1) return;

        std::vector<diagonal_matrix<T> *> split;
        int tailm = __a_mod(a.size, AMBIENT_IB);
        split.reserve(a.nt);
        for(int i = 1; i < a.nt; i++) 
            split.push_back(new diagonal_matrix<T>(AMBIENT_IB));
        split.push_back(new diagonal_matrix<T>(tailm));

        const diagonal_matrix<T>* src = a.data[0];

        if(!a[0].get_data().impl->weak())
        for(int i = 0; i < a.nt; i++){
            diagonal_matrix<T>& dst = *split[i];
            copy_block(src->get_data(), i*AMBIENT_IB, 0, dst.get_data(), 0, 0, dst.num_rows(), 1);
        }

        delete src;
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), split);
    }
    // }}}

    template<class MatrixA, class MatrixB>
    inline void copy_block(const tiles<MatrixB>& in, size_t ii, size_t ij,
                           tiles<MatrixA>& out, size_t oi, size_t oj, 
                           size_t m, size_t n)
    {
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_block(in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
                   out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
                   row.step, col.step);
    }

    template<class Matrix>
    inline void copy_block_s(const tiles<Matrix>& in, size_t ii, size_t ij,
                             tiles<Matrix>& out, size_t oi, size_t oj, 
                             const tiles<Matrix>& alfa, size_t ai, size_t aj,
                             size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= AMBIENT_IB; aj %= AMBIENT_IB;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_block_s(in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
                     out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
                     factor, ai, aj, row.step, col.step);
    }

    template<class Matrix>
    inline void copy_block_sa(const tiles<Matrix>& in, size_t ii, size_t ij,
                              tiles<Matrix>& out, size_t oi, size_t oj, 
                              const tiles<Matrix>& alfa, size_t ai, size_t aj,
                              size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= AMBIENT_IB; aj %= AMBIENT_IB;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_block_sa(in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
                      out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
                      factor, ai, aj, row.step, col.step);
    }

    template<class MatrixA, class MatrixB, class MatrixC>
    inline void gemm(const tiles<MatrixA>& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                std::vector<MatrixC*> ctree; ctree.reserve(a.nt);
                ctree.push_back(&c.tile(i,j));
                size_t rows = c.tile(i,j).num_rows();
                size_t cols = c.tile(i,j).num_cols();
                for(int k = 1; k < a.nt; k++) 
                    ctree.push_back(new MatrixC(rows, cols));
                for(int k = 0; k < a.nt; k++){
                    const MatrixA& ab = a.tile(i, k);
                    const MatrixB& bb = b.tile(k, j);
                    gemm(ab, bb, *ctree[k]);
                }
                __a_reduce(ctree);
                for(int k = 1; k < a.nt; k++) 
                    delete ctree[k];
            }
        }
    }

    template<class MatrixA, class MatrixC, typename T>
    inline void gemm(const tiles<MatrixA>& a, const tiles<diagonal_matrix<T> >& b, tiles<MatrixC>& c){
        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                gemm(a.tile(i,j), b[j], c.tile(i,j));
            }
        }
    }

    template<class MatrixB, class MatrixC, typename T>
    inline void gemm(const tiles<diagonal_matrix<T> >& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                gemm(a[i], b.tile(i,j), c.tile(i,j));
            }
        }
    }

    template<int alfa, int beta, class MatrixA, class MatrixB, class MatrixC>
    inline void gemv(const MatrixA&  a, int ai, int aj, 
                     const MatrixB&  b, int bi, int bj, 
                           MatrixC&  c, int ci, int cj,
                           int m, int n)
    {
        if(m == 0 || n == 0) return;

        for(cross_iterator row(ai,ci,m); !row.end(); ++row){
            std::vector<MatrixC*> ctree;
            for(cross_iterator col(aj,bi,n); !col.end(); ++col){
                MatrixC* part = new MatrixC(row.step, 1);
                gemv<alfa,0>(a.locate(row.first,col.first), a.addr(row.first, col.first),
                             b.locate(col.second,bj), b.addr(col.second, bj),
                             *part, 0,
                             row.step, col.step);
                ctree.push_back(part);
            }
            __a_reduce(ctree);
            add_vectors<beta>(c.locate(row.second, cj), c.addr(row.second, cj), *ctree[0], 0, row.step);
            for(int k = 0; k < ctree.size(); k++) delete ctree[k];
        }
    }
                
    template<class Matrix>
    void orgqr(const tiles<Matrix>&& a, tiles<Matrix>&& q, const tiles<Matrix>&& t){
        int k, m, n;

        for(k = std::min(a.mt, a.nt)-1; k >= 0; k--){
            for(m = q.mt - 1; m > k; m--)
                for(n = 0; n < q.nt; n++)
                    tsmqr<PlasmaNoTrans>(a.tile(m, k).num_cols(), q.tile(k, n), q.tile(m, n), a.tile(m, k), t.tile(m, k));
            
            for(n = 0; n < q.nt; n++)
                ormqr<PlasmaNoTrans>(std::min(a.tile(k,k).num_rows(), a.tile(k,k).num_cols()), a.tile(k, k), t.tile(k, k), q.tile(k, n));
        }
    }

    template<class Matrix>
    void orglq(const tiles<Matrix>&& a, tiles<Matrix>&& q, const tiles<Matrix>&& t){
        int k, m, n;

        for(k = std::min(a.mt, a.nt)-1; k >= 0; k--){
            for(n = q.nt-1; n > k; n--)
                for(m = 0; m < q.mt; m++)
                    tsmlq<PlasmaNoTrans>(a.tile(k, n).num_rows(), q.tile(m, k), q.tile(m, n), a.tile(k, n), t.tile(k, n));
        
            for(m = 0; m < q.mt; m++)
                ormlq<PlasmaNoTrans>(std::min(a.tile(k,k).num_rows(), a.tile(k,k).num_cols()), a.tile(k, k), t.tile(k, k), q.tile(m, k));
        }
    }

    template<size_t OFF = 0, class Matrix>
    void laset2lower(tiles<Matrix>&& a){
        for(size_t j = 0; j < std::min(a.mt, a.nt); j++){
            laset2<PlasmaLower, OFF>(a.tile(j,j));
        
            for(size_t i = j+1; i < a.mt; i++)
                laset2<PlasmaUpperLower, 0>(a.tile(i,j));
        }
    }

    template<size_t OFF = 0, class Matrix>
    void laset2upper(tiles<Matrix>&& a){
        for(size_t j = 1; j < a.nt; j++)
            for(size_t i = 0; i < std::min(j, a.mt); i++)
                laset2<PlasmaUpperLower, 0>(a.tile(i,j));
        
        for(size_t j = 0; j < std::min(a.mt, a.nt); j++)
            laset2<PlasmaUpper, OFF>(a.tile(j,j));
    }

    template<PLASMA_enum LR, class Matrix>
    void orgbr(const tiles<Matrix>& a, tiles<Matrix>& q, const tiles<Matrix>& t){
        if(LR == PlasmaLeft){
            if(num_rows(a) >= num_cols(a)){
                orgqr(a, q, t);
            }else{
                orgqr(a.subset(1, 0, a.mt-1, a.nt), 
                      q.subset(1, 1, q.mt-1, q.nt-1),
                      t.subset(1, 0, t.mt-1, t.nt));
                /*
                   Shift the vectors which define the elementary reflectors one
                   column to the right, and set the first row and column of Q
                   to those of the unit matrix
                   
                    ----------------             ----------------
                   |   |   |   |   |            | 1 | 0 | 0 | 0 |
                    ----------------             ----------------
                   | * |   |   |   |   -- >     | 0 | * |   |   |
                    ----------------             ----------------
                   |   | * |   |   |            | 0 |   | * |   |
                   ----------------             ----------------

                 */
            }
        }else{
            if(num_rows(a) < num_cols(a)){
                orglq(a, q, t);
            }else{
                orglq(a.subset(0, 1, a.mt,   a.nt-1), 
                      q.subset(1, 1, q.mt-1, q.nt-1),
                      t.subset(0, 1, t.mt,   t.nt-1));
                /*
                   Shift the vectors which define the elementary reflectors one
                   row downward, and set the first row and column of P' to
                   those of the unit matrix

                    ----------------             ----------------
                   |   | * |   |   |            | 1 | 0 | 0 | 0 |
                    ----------------             ----------------
                   |   |   | * |   |   -- >     | 0 | * |   |   |
                    ----------------             ----------------
                   |   |   |   | * |            | 0 |   | * |   |
                    ----------------             ----------------
                                                | 0 |   |   | * |
                                                 ----------------
                 */
            }
        }
    }

    template<class Matrix>
    inline void compress_band(tiles<Matrix>& a){
        Matrix* c;
        if(a.num_rows() >= a.num_cols()){
            c = new Matrix(std::min((size_t)(AMBIENT_IB+1),a.num_cols()),a.num_cols());
            copy_band<PlasmaUpper>(a.tile(0,0), *c, 0);
            for(int j = 1; j < a.nt; j++){
                copy_band<PlasmaLower>(a.tile(j-1,j), *c, AMBIENT_IB*j);
                copy_band<PlasmaUpper>(a.tile(j,j),   *c, AMBIENT_IB*j);
            }
        }else{
            c = new Matrix(std::min((size_t)(AMBIENT_IB+1),a.num_rows()),a.num_rows());
            for(int j = 0; j < a.mt-1; j++){
                copy_band<PlasmaLower>(a.tile(j,j),   *c, AMBIENT_IB*j);
                copy_band<PlasmaUpper>(a.tile(j+1,j), *c, AMBIENT_IB*j);
            }
            copy_band<PlasmaLower>(a.tile(a.mt-1,a.mt-1), *c, AMBIENT_IB*(a.mt-1));
        }
        tiles<Matrix> t(c);
        a.swap(t);
    }

    template<class Matrix>
    inline void band(tiles<Matrix> a, tiles<Matrix>& u, tiles<Matrix>& b, tiles<Matrix>& v){
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        resize(u, m, m);
        resize(v, n, n);

        tiles<Matrix> t(a.mt*AMBIENT_IB, a.nt*AMBIENT_IB);
        for(int i = 0; i < std::min(u.mt, u.nt); i++) fill_identity(u.tile(i,i));
        for(int i = 0; i < std::min(v.mt, v.nt); i++) fill_identity(v.tile(i,i));

        if(a.num_rows() >= a.num_cols()){ // upper band diagonal
            for(int k = 0; k < a.nt; k++){
                qr(a.subset(k, k, a.mt-k, 1), t.subset(k, k, t.mt-k, 1));  
            
                ormqr(a.subset(k, k,   a.mt-k, 1),
                      a.subset(k, k+1, a.mt-k, a.nt-k-1), 
                      t.subset(k, k,   t.mt-k, 1));

                if(k+1 == a.nt) break;

                lq(a.subset(k, k+1, 1, a.nt-k-1), 
                   t.subset(k, k+1, 1, t.nt-k-1));
        
                ormlq(a.subset(k,   k+1, 1,        a.nt-k-1),
                      a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1),
                      t.subset(k,   k+1, 1,        t.nt-k-1));
            }
            orgqr(a, u, t);
            orglq(a.subset(0, 1, a.mt,   a.nt-1), 
                  v.subset(1, 1, v.mt-1, v.nt-1),
                  t.subset(0, 1, t.mt,   t.nt-1));

            laset2upper(a.subset(0, 1, a.mt, a.nt-1));
            laset2lower(a);
        }else{ // lower band diagonal
            for(int k = 0; k < a.mt; k++){
                lq(a.subset(k, k, 1, a.nt-k), 
                   t.subset(k, k, 1, t.nt-k));
            
                ormlq(a.subset(k,   k, 1,        a.nt-k),
                      a.subset(k+1, k, a.mt-k-1, a.nt-k),
                      t.subset(k,   k, 1,        t.nt-k));

                if(k+1 == a.mt) break;

                qr(a.subset(k+1, k, a.mt-k-1, 1),
                   t.subset(k+1, k, t.mt-k-1, 1));
        
                ormqr(a.subset(k+1, k,   a.mt-k-1, 1),
                      a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1),
                      t.subset(k+1, k,   t.mt-k-1, 1));
            }
            orglq(a, v, t);
            orgqr(a.subset(1, 0, a.mt-1, a.nt), 
                  u.subset(1, 1, u.mt-1, u.nt-1),
                  t.subset(1, 0, t.mt-1, t.nt));
            laset2lower(a.subset(1, 0, a.mt-1, a.nt));
            laset2upper(a);
        }

        b.swap(a);
    }

    template<class Matrix>
    inline void shuffle(tiles<Matrix>& a){
        size_t size = a.mt*a.nt;
        for(int i = 0; i < size; i++){
            int dest = (int)(drand48()*(double)size);
            std::swap(a.data[i], a.data[dest]);
        }
    }

    template<class SMatrix, class DiagonalMatrix, class Matrix>
    void plabrd(tiles<SMatrix>&& a, Matrix& say, Matrix& sax, DiagonalMatrix& d, DiagonalMatrix& e, 
                 DiagonalMatrix& tq, DiagonalMatrix& tp, 
                 tiles<Matrix>& x, tiles<Matrix>& y)
    {
        if(num_rows(a) >= num_cols(a)){
       
            int m = num_rows(a);
            int n = num_rows(a);

            merge(x); Matrix& sx = x[0];
            merge(y); Matrix& sy = y[0];
        
            for(int i = 0; i < AMBIENT_IB; ++i)
            {
                labrd_update_col(say, sax, sy, sx, tq, d, i);
                gemv<1,0>(transpose(a), AMBIENT_IB, i, say, i, i, sy, AMBIENT_IB, i, n-AMBIENT_IB, m-i);
                labrd_reduce_col(say, sax, sy, sx, i);
                scale(sy, i+1, i, tq, i, i);
            
                labrd_update_row(say, sax, sy, sx, tp, e, i);
                labrd_reduce_row(say, sax, sy, sx, i);
                gemv<1,1>(a, AMBIENT_IB, i+1, transpose(sax), i+1, i, sx, AMBIENT_IB, i, m-AMBIENT_IB, n-i-1);
                scale(sx, i+1, i, tp, i, i);
            }
            /* {{{ explicit implementation
            for(int i = 0; i < AMBIENT_IB; ++i){
        
                int ri  = m-i;   //std::min(m-i, i*nb);
                int rj  = n-i-1; //std::min(n-i-1, (i+1)*nb);
                int rij = m-i-1; //std::min(m-i-1, (i+1)*nb);
        
                gemv<-1,1>(say,                              i, 0, //
                           transpose(y.subset(0, 0, 1, 1)),  0, i, // only 2 vertical blocks
                           say,                              i, i, //
                           ri, i);                                 // can be grouped
                gemv<-1,1>(x,                                i, 0, //
                           say,                              0, i, // only 1 vertical block
                           say,                              i, i, //
                           ri, i);                                 //
                                                                   //
                larfg<PlasmaNoTrans>(say, tq, d, i);               //
                
                // --------------------- BIG ONE -------------------------
                gemv<1,0>(transpose(a),                            i+1, i,
                          say,                                     i,   i,
                          y,                                       i+1, i, // 0
                          rj, ri);
                // -------------------------------------------------------
                gemv<1,0>(transpose(say),                          0,   i, // {
                          say,                                     i,   i, //
                          y.subset(0, 0, 1, 1),                    0,   i, //
                          i, ri);                                          // can be groupped
                gemv<-1,1>(y,                                      i+1, 0, //
                           y.subset(0, 0, 1, 1),                   0,   i, //
                           y,                                      i+1, i, // } 0
                           rj, i);
                gemv<1,0>(transpose(x),                            0,   i, // {
                          say,                                     i,   i, //
                          y.subset(0, 0, 1, 1),                    0,   i, //
                          i, ri);                                          // can be groupped
                gemv<-1,1>(transpose(sax),                         i+1, 0, // 
                           y.subset(0, 0, 1, 1),                   0,   i, //
                           y,                                      i+1, i, // } 0
                           rj, i);
        
                scale(y, i+1, i, tq, i, i);

                // synchronizing stripes
                copy_block(say, 0, 0, sax, 0, 0, num_rows(sax), num_cols(say));
        
                gemv<-1,1>(y,                                      i+1, 0, //
                           transpose(sax),                         0,   i, // only 2 horizontal blocks
                           transpose(sax),                         i+1, i, // can be grouped
                           rj, i+1);                                       //
                gemv<-1,1>(transpose(sax),                         i+1, 0, //
                           transpose(x.subset(0, 0, 1, 1)),        0,   i, //
                           transpose(sax),                         i+1, i, //
                           rj, i);                                         //
                                                                           //
                larfg<PlasmaTrans>(sax, tp, e, i);                         //
        
                copy_block(sax, 0, 0, say, 0, 0, num_rows(sax), num_cols(say));
                // synchronizing stripes

                // --------------------- BIG ONE -------------------------
                gemv<1,0>(a,                                       i+1, i+1,
                          transpose(sax),                          i+1, i,
                          x,                                       i+1, i,   // 1
                          rij, rj);                     
                // -------------------------------------------------------
                gemv<1,0>(transpose(y),                            0,   i+1, // {
                          transpose(sax),                          i+1, i,   //
                          x.subset(0, 0, 1, 1),                    0,   i,   //
                          i+1, rj);                                          // can be groupped
                gemv<-1,1>(say,                                    i+1, 0,   //
                           x.subset(0, 0, 1, 1),                   0,   i,   //
                           x,                                      i+1, i,   // } 1
                           rij, i+1);
                gemv<1,0>(sax,                                     0,   i+1, // {
                          transpose(sax),                          i+1, i,   //
                          x.subset(0, 0, 1, 1),                    0,   i,   //
                          i, rj);                                            // can be groupped
                gemv<-1,1>(x,                                      i+1, 0,   //
                           x.subset(0, 0, 1, 1),                   0,   i,   //
                           x,                                      i+1, i,   // } 1
                           rij, i);
        
                scale(x, i+1, i, tp, i, i);
            }
            }}} */

            split_d(x);
            split_d(y);

        }else{}
    }

    template<class Matrix, class DiagonalMatrix>
    void pgebd2(Matrix& a, DiagonalMatrix& d, DiagonalMatrix& e, DiagonalMatrix& tq, DiagonalMatrix& tp){
        gebd2(a, d, e, tq, tp);
    }

    template<class Matrix, class DiagonalMatrix>
    void pgebrd(tiles<Matrix>& a, tiles<DiagonalMatrix>& d, tiles<DiagonalMatrix>& e, tiles<Matrix>& u, tiles<Matrix>& v)
    {
        int m = num_rows(a);
        int n = num_cols(a);
        int k = std::min(m,n);
        
        tiles<DiagonalMatrix> tp(k);
        tiles<DiagonalMatrix> tq(k);
        tiles<Matrix> x(m, AMBIENT_IB);
        tiles<Matrix> y(n, AMBIENT_IB);
        resize(d, k, k);
        resize(e, k, k);
        
        for(k = 0; k < std::min(a.mt, a.nt)-1; ++k)
        {
            Matrix& say = *merge(a.subset(k, k, a.mt-k, 1));
            Matrix& sax = *merge(a.subset(k, k, 1, a.nt-k));

            plabrd(a.subset(k, k, a.mt-k, a.nt-k), say, sax, d[k], e[k], tq[k], tp[k], x, y);

            for(int i = k+1; i < a.mt; i++)
                copy_block(say, (i-k)*AMBIENT_IB, 0, a.tile(i,k), 0, 0, a.tile(i,k).num_rows(), a.tile(i,k).num_cols());
            for(int j = k+1; j < a.nt; j++)
                copy_block(sax, 0, (j-k)*AMBIENT_IB, a.tile(k,j), 0, 0, a.tile(k,j).num_rows(), a.tile(k,j).num_cols());

            delete &say;
            delete &sax;

            a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1) -= a.subset(k+1, k, a.mt-k-1, 1) * transpose(y.subset(1, 0, y.mt-k-1, 1));
            a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1) -= x.subset(  1, 0, x.mt-k-1, 1) * a.subset(k, k+1, 1, a.nt-k-1);
        }
        Matrix* tail = merge(a.subset(k, k, a.mt-k, a.nt-k));
        pgebd2(*tail, d[k], e[k], tq[k], tp[k]);
        delete tail;
    }

    template<class Matrix, class DiagonalMatrix>
    inline void svd_mod(const tiles<Matrix>& a, tiles<Matrix>& u, tiles<Matrix>& v, tiles<DiagonalMatrix>& s){
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        size_t k = std::min(m,n);
        tiles<DiagonalMatrix> e(k);

        s.resize(k, k);
        u.resize(m, m); 
        tiles<Matrix> u1;
        tiles<Matrix> u2(new Matrix(m,m));
        tiles<Matrix> s1;
        tiles<Matrix> v2(new Matrix(n,n));
        tiles<Matrix> v1;
        v.resize(n, n); 

        band(a, u1, s1, v1);

#ifdef GBBRD
        compress_band(s1);
        merge(s1); merge(s); merge(e);
        gbbrd(s1[0], s[0], e[0], u2[0], v2[0]);
#elif defined(GEBRD)
        merge(s1); merge(s); merge(e);
        gebrd(s1[0], s[0], e[0], u2[0], v2[0]);
#else
        pgebrd(s1, s, e, u2, v2);
        merge(s); merge(e);
#endif
        bdsqr(s[0], e[0], u2[0], v2[0]);
        split_d(s); split_d(u2); split_d(v2);

        ambient::numeric::gemm(u1, u2, u);
        ambient::numeric::gemm(v2, v1, v);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void svd(tiles<Matrix> a, tiles<Matrix>& u, tiles<Matrix>& v, tiles<DiagonalMatrix>& s){
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        size_t k = std::min(m,n);
        resize(u, m, k);
        resize(s, k, k);
        resize(v, k, n);

        merge(a); merge(u); merge(v); merge(s);
        svd(a[0], u[0], v[0], s[0]);
        split_d(u); split_d(v); split_d(s);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void heev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        merge(a); merge(evecs); merge(evals);
        heev(a[0], evecs[0], evals[0]);
        split_d(a); split_d(evecs); split_d(evals);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void syev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        heev(a, evecs, evals);
    }

    template<class Matrix>
    void ormqr(tiles<Matrix>&& a, tiles<Matrix>&& b, tiles<Matrix>&& t){
        int k, m, n;
    
       // PlasmaLeft / PlasmaTrans
        for(k = 0; k < std::min(a.mt, a.nt); k++){
            size_t kmin = std::min(a.tile(k,k).num_rows(), a.tile(k,k).num_cols());
            for(n = 0; n < b.nt; n++){
                ormqr<PlasmaTrans>(kmin, a.tile(k, k), t.tile(k, k), b.tile(k, n));
            }
            for(m = k+1; m < b.mt; m++){
                for(n = 0; n < b.nt; n++)
                    tsmqr<PlasmaTrans>(kmin, b.tile(k, n), b.tile(m, n), a.tile(m, k), t.tile(m, k));
            }
        }
    }

    template<class Matrix>
    void ormlq(tiles<Matrix>&& a, tiles<Matrix>&& b, tiles<Matrix>&& t){
        int k, m, n;
        
        // PlasmaRight / PlasmaTrans
        for(k = 0; k < std::min(a.mt, a.nt); k++){
            size_t kmin = std::min(a.tile(k,k).num_rows(), a.tile(k,k).num_cols());
            for(m = 0; m < b.mt; m++){
                ormlq<PlasmaTrans>(kmin, a.tile(k, k), t.tile(k, k), b.tile(m, k));
            }
            for(n = k+1; n < b.nt; n++){
                for(m = 0; m < b.mt; m++){
                    tsmlq<PlasmaTrans>(kmin, b.tile(m, k), b.tile(m, n), a.tile(k, n), t.tile(k, n));
                }
            }
        }
    }

    template<class Matrix>
    inline void qr(tiles<Matrix>&& a, tiles<Matrix>&& t){
        int k, m, n;
        
        for(k = 0; k < std::min(a.mt, a.nt); k++){
            geqrt(a.tile(k, k), t.tile(k, k));
        
            for(n = k+1; n < a.nt; n++)
                ormqr<PlasmaTrans>(a.tile(k,n).num_rows(), a.tile(k, k), t.tile(k, k), a.tile(k, n));
            
            for(m = k+1; m < a.mt; m++){
                tsqrt(a.tile(k, k), a.tile(m, k), t.tile(m, k));
        
                for(n = k+1; n < a.nt; n++)
                    tsmqr<PlasmaTrans>(AMBIENT_IB, a.tile(k, n), a.tile(m, n), a.tile(m, k), t.tile(m, k));
            }
        }
    }


    template<class Matrix>
    inline void qr(tiles<Matrix> a, tiles<Matrix>& q, tiles<Matrix>& r){
        int am = num_rows(a);
        int an = num_cols(a);
        int ak = std::min(am,an);
        resize(q, am, an);
        resize(r, ak, an);
        
        tiles<Matrix> t(a.mt*AMBIENT_IB, a.nt*AMBIENT_IB);
        qr(a, t);
        
        // restoring R from A //
        for(int j = 0; j < a.nt; j++)
        for(int i = 0; i < j && i < std::min(a.mt, a.nt); i++)
            copy(a.tile(i,j), r.tile(i,j));
        
        for(int k = 0; k < std::min(a.mt, a.nt); k++)
            copy_rt(a.tile(k,k), r.tile(k,k));
        
        // restoring Q from T //
        // do we need to memset Q first?
        for(int i = 0; i < std::min(a.mt, a.nt); i++)
            fill_identity(q.tile(i,i));
       
        orgqr(a, q, t); 
        resize(q, am, ak);
    }

    template<class Matrix>
    inline void lq(tiles<Matrix>&& a, tiles<Matrix>&& t){
        int k, m, n;
        
        for(k = 0; k < std::min(a.mt, a.nt); k++) {
            gelqt(a.tile(k, k), t.tile(k, k));
        
            for(m = k+1; m < a.mt; m++)
                ormlq<PlasmaTrans>(a.tile(m, k).num_cols(), a.tile(k, k), t.tile(k, k), a.tile(m, k));
        
            for(n = k+1; n < a.nt; n++){
                tslqt(a.tile(k, k), a.tile(k, n), t.tile(k, n));
        
                for(m = k+1; m < a.mt; m++)
                    tsmlq<PlasmaTrans>(AMBIENT_IB, a.tile(m, k), a.tile(m, n), a.tile(k, n), t.tile(k, n));
            }
        }
    }

    template<class Matrix>
    inline void lq(tiles<Matrix> a, tiles<Matrix>& l, tiles<Matrix>& q){
        int am = num_rows(a);
        int an = num_cols(a);
        int ak = std::min(am,an);
        resize(l, am, ak);
        resize(q, am, an); // instead of ak
        
        tiles<Matrix> t(a.mt*AMBIENT_IB, a.nt*AMBIENT_IB);
        lq(a, t);
        
        // restoring L from A //
        for(int j = 0; j < std::min(a.mt, a.nt); ++j)
        for(int i = j+1; i < a.mt; ++i)
            copy(a.tile(i,j), l.tile(i,j));
        
        for(int k = 0; k < std::min(a.mt, a.nt); k++)
            copy_lt(a.tile(k,k), l.tile(k,k));
        
        // restoring Q from T //
        // do we need to memset Q first?
        for(int i = 0; i < std::min(a.mt, a.nt); i++)
            fill_identity(q.tile(i,i));
       
        orglq(a, q, t); 
        resize(q, ak, an);
    }

    template<class Matrix> 
    inline void resize(tiles<Matrix>& a, size_t m, size_t n){ 
        if(a.num_rows() == m && a.num_cols() == n) return;
        tiles<Matrix> r(m, n);

        int mb_min = std::min(r.mt, a.mt);
        int nb_min = std::min(r.nt, a.nt);
        
        for(int j = 0; j < nb_min; j++){
            for(int i = 0; i < mb_min; i++){
                r.tile(i,j) = a.tile(i,j);
            }
        }
        size_t margin = AMBIENT_IB;
        if(r.mt <= a.mt) margin = __a_mod(m, AMBIENT_IB);
        for(int j = 0; j < nb_min; j++){
            Matrix& block = r.tile(mb_min-1,j);
            resize(block, margin, block.num_cols());
        }
        margin = AMBIENT_IB;
        if(r.nt <= a.nt) margin = __a_mod(n, AMBIENT_IB);
        for(int i = 0; i < mb_min; i++){
            Matrix& block = r.tile(i,nb_min-1);
            resize(block, block.num_rows(), margin);
        }
        swap(a, r);
    }

    template<typename T> 
    inline void resize(tiles<diagonal_matrix<T> >& a, size_t m, size_t n){
        if(a.num_rows() == m) return;
        tiles<diagonal_matrix<T> > r(m);

        int nb_min = std::min(r.nt, a.nt);
        
        for(int i = 0; i < nb_min; i++){
            r[i] = a[i];
        }
        if(r.nt > a.nt){
            resize(r[nb_min-1], AMBIENT_IB, AMBIENT_IB);
        }else{
            size_t margin = __a_mod(m, AMBIENT_IB);
            resize(r[r.nt-1], margin, margin);
        }
        swap(a, r);
    }

    template<typename T> 
    inline void sqrt_inplace(tiles<diagonal_matrix<T> >& a){
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            sqrt_inplace(a[i]);
        }
    }

    template<typename T>
    inline void exp_inplace(tiles<diagonal_matrix<T> >& a, const T& alfa = 1.){
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            exp_inplace(a[i], alfa);
        }
    }

    template<class Matrix>
    inline tiles<Matrix> exp(const tiles<Matrix>& a, const value_type& alfa = 1.){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
    }

    template<class Matrix> 
    inline scalar_type trace(const tiles<Matrix>& a){
        int size = std::min(a.nt, a.mt);
        scalar_type result(0.);
        std::vector<scalar_type> parts;
        parts.reserve(size);
        
        for(int k = 0; k < size; k++) parts.push_back(trace(a.tile(k, k)));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline real_type norm_square(const tiles<Matrix>& a){
        int size = a.nt*a.mt;
        real_type result(0.);
        std::vector<real_type> parts;
        parts.reserve(size);
        
        for(int k = 0; k < size; k++) parts.push_back(norm_square(a[k]));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline scalar_type overlap(const tiles<Matrix>& a, const tiles<Matrix>& b){
        int size = a.nt*a.mt;
        scalar_type result(0.);
        std::vector<scalar_type> parts;
        parts.reserve(size);
        
        for(int k = 0; k < size; k++) parts.push_back(overlap(a[k], b[k]));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }
        
    template<class Matrix>
    inline void swap(tiles<Matrix>& x, tiles<Matrix>& y){ 
        x.swap(y);                     
    }

    template<class Matrix>
    inline const tiles<Matrix>& conj(const tiles<Matrix>& a){
        //a.conj();
        return a;
    }

    template<class Matrix>
    inline void conj_inplace(tiles<Matrix>& a){
        // gs (doubles)
        // does nothing for now
    }

    template<class Matrix>
    inline void transpose_inplace(tiles<Matrix>& a){
        std::vector<Matrix*> t;
        t.reserve(a.mt*a.nt);
        for(int i = 0; i < a.mt; i++){
            for(int j = 0; j < a.nt; j++){
                Matrix& block = a.tile(i,j);
                transpose_inplace(block);
                t.push_back(&block);
            }
        }
        std::swap(a.data, t);
        std::swap(a.rows, a.cols);
        std::swap(a.mt,   a.nt);
    }

    template<class Matrix>
    inline tiles<transpose_view<Matrix> > transpose(const tiles<Matrix>& a){
        tiles<transpose_view<Matrix> > t;
        std::vector<transpose_view<Matrix>*> data;
        data.reserve(a.mt*a.nt);
        for(int i = 0; i < a.mt; i++)
            for(int j = 0; j < a.nt; j++)
                data.push_back(new transpose_view<Matrix>(a.tile(i,j)));
        std::swap(t.data, data);
        t.cols = a.num_rows();
        t.rows = a.num_cols();
        t.nt   = a.mt;
        t.mt   = a.nt;
        return t;
    }

    template<class Matrix>
    inline void adjoint_inplace(tiles<Matrix>& a){
        transpose_inplace(a);
    }

    template<class Matrix>
    inline tiles<transpose_view<Matrix> > adjoint(const tiles<Matrix>& a){
        return transpose(a);
    }

    template<class Matrix, class G>
    inline void generate(tiles<Matrix>& a, G g){
        generate(a);
    }

    template<class Matrix>
    inline void generate(tiles<Matrix>& a){
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            fill_random(a[i]);
        }
    }

    template<class Matrix>
    inline void remove_rows(tiles<Matrix>& a, size_type i, difference_type k){
        assert(false); printf("ERROR: NOT TESTED (BLOCKED REMOVE ROWS)\n");
    }

    template<class Matrix>
    inline void remove_cols(tiles<Matrix>& a, size_type j, difference_type k){
        assert(false); printf("ERROR: NOT TESTED (BLOCKED REMOVE COLS)\n");
    }

    template <class MatrixA, class MatrixB>
    inline void add_inplace(tiles<MatrixA>& lhs, const tiles<MatrixB>& rhs){
        int size = lhs.data.size();
        for(int i = 0; i < size; i++){
            lhs[i] += rhs[i];
        }
    }

    template <class MatrixA, class MatrixB>
    inline void sub_inplace(tiles<MatrixA>& lhs, const tiles<MatrixB>& rhs){ 
        int size = lhs.data.size();
        for(int i = 0; i < size; i++){
            lhs[i] -= rhs[i];
        }
    }

    template <class MatrixA, class MatrixB>
    inline void mul_inplace(tiles<MatrixA>& a, const tiles<MatrixB>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM BLOCKED INPLACE)\n");
    }

    template <class Matrix>
    inline void mul_inplace(tiles<Matrix>& a, const scalar_type& rhs) { 
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            a[i] *= rhs;
        }
    }

    template <class Matrix>
    inline void div_inplace(tiles<Matrix>& a, const scalar_type& rhs){
        int size = a.data.size();
        for(int i = 0; i < size; i++)
            a[i] /= rhs;
    }
 
    template <class Matrix>
    inline void save(const tiles<Matrix>& a, size_t tag){
        split_d(a);
        int size = a.data.size();
        for(int i = 0; i < size; ++i)
           save(a[i], (tag+i));//tag is done over blocks      
    }

    template <class Matrix>
    inline void load(tiles<Matrix>& a, size_t tag){
        split_d(a);
        int size = a.data.size();
        for(int i = 0; i < size; ++i)
           load(a[i], (tag+i));        
    }

    template <class Matrix>
    std::ostream& operator << (std::ostream& o, const tiles<Matrix>& a){
        int size = a.data.size();
        for(int i = 0; i < size; i++)
            o << a[i];
        return o;
    }

    template <class Matrix> 
    inline tiles<Matrix> operator + (tiles<Matrix> lhs, const tiles<Matrix>& rhs){ 
        return (lhs += rhs);
    }

    template <class Matrix> 
    inline tiles<Matrix> operator - (tiles<Matrix> lhs, const tiles<Matrix>& rhs){ 
        return (lhs -= rhs);
    }

    template <class Matrix, class MatrixB>
    inline tiles<matrix<value_type> > operator * (const tiles<Matrix>& lhs, const tiles<MatrixB>& rhs){
        tiles<matrix<value_type> > ret(lhs.num_rows(), rhs.num_cols());
        gemm(lhs, rhs, ret);
        return ret; 
    }

    template<class Matrix, typename T> 
    inline const tiles<Matrix> operator * (tiles<Matrix> lhs, const T& rhs){ 
        return (lhs *= rhs); 
    }

    template<class Matrix, typename T> 
    inline const tiles<Matrix> operator * (const T& lhs, tiles<Matrix> rhs){ 
        return (rhs *= lhs); 
    }

    template<class Matrix> inline size_type num_rows(const tiles<Matrix>& a){ 
        return a.num_rows();
    }

    template<class Matrix> inline size_type num_cols(const tiles<Matrix>& a){
        return a.num_cols();
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type 
#endif
