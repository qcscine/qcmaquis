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
    bool is_hermitian(const tiles<Matrix>& m){
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
                copy(dst, 0, 0, *src, i*AMBIENT_IB, j*AMBIENT_IB, dst.num_rows(), dst.num_cols());
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
            copy(dst.get_data(), 0, 0, src->get_data(), i*AMBIENT_IB, 0, dst.num_rows(), 1);
        }

        delete src;
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), split);
    }

    // {{{ normal merge / split
    
    template<class Matrix>
    inline void merge(const tiles<Matrix>& a){
        if(a.data.size() == 1) return;

        std::vector<Matrix*> merge; 
        merge.push_back(new Matrix(a.rows, a.cols));

        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                const Matrix* src = &a.tile(i,j);
                if(!src->impl->weak())
                copy(*merge[0], i*AMBIENT_IB, j*AMBIENT_IB, *src, 0, 0, src->num_rows(), src->num_cols());
                delete src;
            }
        }
        std::swap(const_cast<std::vector<Matrix*>&>(a.data), merge);
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
                copy(dst, 0, 0, *src, i*AMBIENT_IB, j*AMBIENT_IB, dst.num_rows(), dst.num_cols());
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

        std::vector<diagonal_matrix<T> *> merge; 
        merge.push_back(new diagonal_matrix<T>(a.size));

        for(int i = 0; i < a.nt; i++){
            const diagonal_matrix<T>* src = a.data[i];
            if(!src->get_data().impl->weak())
            copy(merge[0]->get_data(), i*AMBIENT_IB, 0, src->get_data(), 0, 0, src->num_rows(), 1);
            delete src;
        }
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), merge);
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
            copy(dst.get_data(), 0, 0, src->get_data(), i*AMBIENT_IB, 0, dst.num_rows(), 1);
        }

        delete src;
        std::swap(const_cast<std::vector<diagonal_matrix<T>*>&>(a.data), split);
    }
    // }}}

    template<class Matrix>
    inline void copy(tiles<Matrix>& out, size_t oi, size_t oj, 
                     const tiles<Matrix>& in, size_t ii, size_t ij,
                     size_t m, size_t n)
    {
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy(out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
             in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
             row.step, col.step);
    }

    template<class Matrix>
    inline void copy_s(tiles<Matrix>& out, size_t oi, size_t oj, 
                       const tiles<Matrix>& in, size_t ii, size_t ij,
                       const tiles<Matrix>& alfa, size_t ai, size_t aj,
                       size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= AMBIENT_IB; aj %= AMBIENT_IB;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_s(out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
               in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
               factor, ai, aj, row.step, col.step);
    }

    template<class Matrix>
    inline void copy_sa(tiles<Matrix>& out, size_t oi, size_t oj, 
                        const tiles<Matrix>& in, size_t ii, size_t ij,
                        const tiles<Matrix>& alfa, size_t ai, size_t aj,
                        size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= AMBIENT_IB; aj %= AMBIENT_IB;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_sa(out.locate(row.first, col.first), row.first%AMBIENT_IB, col.first%AMBIENT_IB, 
                in.locate(row.second, col.second), row.second%AMBIENT_IB, col.second%AMBIENT_IB,
                factor, ai, aj, row.step, col.step);
    }

    template<class MatrixA, class MatrixB, class MatrixC>
    inline void gemm(const tiles<MatrixA>& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        split_d(a); split_d(b); split_d(c);
        
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
        split_d(a); split_d(b); split_d(c);
        
        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                gemm(a.tile(i,j), b[j], c.tile(i,j));
            }
        }
    }

    template<class MatrixB, class MatrixC, typename T>
    inline void gemm(const tiles<diagonal_matrix<T> >& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        split_d(a); split_d(b); split_d(c);
        
        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                gemm(a[i], b.tile(i,j), c.tile(i,j));
            }
        }
    }

    template<PLASMA_enum UL, class Matrix>
    void laset2(tiles<Matrix>&& a, const value_type alfa = 0.0){
        size_t i, j;

        if(UL == PlasmaLower){
           for(j = 0; j < std::min(a.mt, a.nt); j++){
               laset2<PlasmaLower>(a.tile(j,j), alfa);
        
               for(i = j+1; i < a.mt; i++)
                   laset2<PlasmaUpperLower>(a.tile(i,j), alfa);
           }
        }else{
           for(j = 1; j < a.nt; j++)
               for(i = 0; i < std::min(j, a.mt); i++)
                   laset2<PlasmaUpperLower>(a.tile(i,j), alfa);
           
           for(j = 0; j < std::min(a.mt, a.nt); j++)
               laset2<PlasmaUpper>(a.tile(j,j), alfa);
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
    void gebrd_ge2tb(tiles<Matrix>& a, tiles<Matrix>& t){ // band reduction
        int k;
        
        if(a.num_rows() >= a.num_cols()){
            for(k = 0; k < a.nt; k++){
                qr(a.subset(k, k, a.mt-k, 1), t.subset(k, k, t.mt-k, 1));  
            
                ormqr(a.subset(k, k,   a.mt-k, 1),
                      a.subset(k, k+1, a.mt-k, t.nt-k-1), 
                      t.subset(k, k,   t.mt-k, 1));
                if(k+1 < a.nt){
                   lq(a.subset(k, k+1, 1, a.nt-k-1), 
                      t.subset(k, k+1, 1, t.nt-k-1));
        
                   ormlq(a.subset(k,   k+1, 1,        a.nt-k-1),
                         a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1),
                         t.subset(k,   k+1, 1,        t.nt-k-1));
                }
            }
        }else{
            for(k = 0; k < a.mt; k++){
                lq(a.subset(k, k, 1, a.nt-k), 
                   t.subset(k, k, 1, t.nt-k));
            
                ormlq(a.subset(k,   k, 1,        a.nt-k),
                      a.subset(k+1, k, a.mt-k-1, a.nt-k),
                      t.subset(k,   k, 1,        t.nt-k));
                if(k+1 < a.mt){
                   qr(a.subset(k+1, k, a.mt-k-1, 1),
                      t.subset(k+1, k, t.mt-k-1, 1));
        
                   ormqr(a.subset(k+1, k,   a.mt-k-1, 1),
                         a.subset(k+1, k+1, a.mt-k-1, a.nt-k-1),
                         t.subset(k+1, k,   t.mt-k-1, 1));
                }
            }
        }
    }

    template<class Matrix>
    inline void gebrd(tiles<Matrix> a, tiles<Matrix>& u, tiles<Matrix>& b, tiles<Matrix>& vt){
        split_d(a);

        size_t m = num_rows(a);
        size_t n = num_cols(a);
        size_t k = std::min(m,n);
        resize(u,  m, n);
        resize(vt, m, n);

        tiles<Matrix> t(a.mt*AMBIENT_IB, a.nt*AMBIENT_IB);
        gebrd_ge2tb(a, t);

        for(int i = 0; i < std::min(u.mt, u.nt); i++) fill_identity(u.tile(i,i));
        for(int i = 0; i < std::min(vt.mt, vt.nt); i++) fill_identity(vt.tile(i,i));

        orgbr<PlasmaLeft>(a, u, t);
        orgbr<PlasmaRight>(a, vt, t);

        if(num_rows(a) < num_cols(a)){
            laset2<PlasmaLower>(a.subset(1, 0, a.mt-1, a.nt));
            laset2<PlasmaUpper>(a);
        }else{
            laset2<PlasmaUpper>(a.subset(0, 1, a.mt, a.nt-1));
            laset2<PlasmaLower>(a);
        }

        resize(u,  m, k);
        resize(a,  k, k);
        resize(vt, k, n);

        b.swap(a);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void svd(tiles<Matrix> a, tiles<Matrix>& u, tiles<Matrix>& vt, tiles<DiagonalMatrix>& s){
        size_t m = num_rows(a);
        size_t n = num_cols(a);
        size_t k = std::min(m,n);
        resize(u,  m, k);
        resize(s,  k, k);
        resize(vt, k, n);

        merge(a); merge(u); merge(vt); merge(s);
        svd(a[0], u[0], vt[0], s[0]);
        split(u); split(vt); split(s);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void heev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        merge(a); merge(evecs); merge(evals);
        heev(a[0], evecs[0], evals[0]);
        split(a); split(evecs); split(evals);
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
        split_d(a); split_d(q); split_d(r);
        
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
            copy(r.tile(i,j), a.tile(i,j));
        
        for(int k = 0; k < std::min(a.mt, a.nt); k++)
            copy_rt(r.tile(k,k), a.tile(k,k));
        
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
        split_d(a); split_d(l); split_d(q);
        
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
            copy(l.tile(i,j), a.tile(i,j));
        
        for(int k = 0; k < std::min(a.mt, a.nt); k++)
            copy_lt(l.tile(k,k), a.tile(k,k));
        
        // restoring Q from T //
        // do we need to memset Q first?
        for(int i = 0; i < std::min(a.mt, a.nt); i++)
            fill_identity(q.tile(i,i));
       
        orglq(a, q, t); 
        resize(q, ak, an);
    }

    template<class Matrix> 
    inline void resize(tiles<Matrix>& m, size_t rows, size_t cols){ 
        if(m.num_rows() == rows && m.num_cols() == cols) return;
        tiles<Matrix> r(rows, cols);

        split_d(m); split_d(r);
        int mb_min = std::min(r.mt, m.mt);
        int nb_min = std::min(r.nt, m.nt);
        
        for(int j = 0; j < nb_min; j++){
            for(int i = 0; i < mb_min; i++){
                r.tile(i,j) = m.tile(i,j);
            }
        }
        size_t margin = AMBIENT_IB;
        if(r.mt <= m.mt) margin = __a_mod(rows, AMBIENT_IB);
        for(int j = 0; j < nb_min; j++){
            Matrix& block = r.tile(mb_min-1,j);
            resize(block, margin, block.num_cols());
        }
        margin = AMBIENT_IB;
        if(r.nt <= m.nt) margin = __a_mod(cols, AMBIENT_IB);
        for(int i = 0; i < mb_min; i++){
            Matrix& block = r.tile(i,nb_min-1);
            resize(block, block.num_rows(), margin);
        }
        swap(m, r);
    }

    template<typename T> 
    inline void resize(tiles<diagonal_matrix<T> >& m, size_t rows, size_t cols){
        if(m.num_rows() == rows) return;
        tiles<diagonal_matrix<T> > r(rows);

        split_d(m); split_d(r);
        int nb_min = std::min(r.nt, m.nt);
        
        for(int i = 0; i < nb_min; i++){
            r[i] = m[i];
        }
        if(r.nt > m.nt){
            resize(r[nb_min-1], AMBIENT_IB, AMBIENT_IB);
        }else{
            size_t margin = __a_mod(rows, AMBIENT_IB);
            resize(r[r.nt-1], margin, margin);
        }
        swap(m, r);
    }

    template<typename T> 
    inline void sqrt_inplace(tiles<diagonal_matrix<T> >& m){
        split_d(m);
        int size = m.data.size();
        for(int i = 0; i < size; i++){
            sqrt_inplace(m[i]);
        }
    }

    template<typename T>
    inline void exp_inplace(tiles<diagonal_matrix<T> >& m, const T& alfa = 1.){
        split_d(m);
        int size = m.data.size();
        for(int i = 0; i < size; i++){
            exp_inplace(m[i], alfa);
        }
    }

    template<class Matrix>
    inline tiles<Matrix> exp(const tiles<Matrix>& m, const value_type& alfa = 1.){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
    }

    template<class Matrix> 
    inline scalar_type trace(const tiles<Matrix>& m){
        split_d(m);
        int size = std::min(m.nt, m.mt);
        scalar_type result(0.);
        std::vector<scalar_type> parts;
        parts.reserve(size);
        
        for(int k = 0; k < size; k++) parts.push_back(trace(m.tile(k, k)));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline real_type norm_square(const tiles<Matrix>& m){ 
        split_d(m);
        int size = m.nt*m.mt;
        real_type result(0.);
        std::vector<real_type> parts;
        parts.reserve(size);
        
        for(int k = 0; k < size; k++) parts.push_back(norm_square(m[k]));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline scalar_type overlap(const tiles<Matrix>& a, const tiles<Matrix>& b){
        split_d(a); split_d(b);
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
    inline const tiles<Matrix>& conj(const tiles<Matrix>& m){
        //m.conj();
        return m;
    }

    template<class Matrix>
    inline void conj_inplace(tiles<Matrix>& m){
        // gs (doubles)
        // does nothing for now
    }

    template<class Matrix>
    inline void transpose_inplace(tiles<Matrix>& a){
        split_d(a);
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
        split_d(a);
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
    inline void adjoint_inplace(tiles<Matrix>& m){
        transpose_inplace(m);
    }

    template<class Matrix>
    inline tiles<transpose_view<Matrix> > adjoint(const tiles<Matrix>& m){
        return transpose(m);
    }

    template<class Matrix, class G>
    inline void generate(tiles<Matrix>& m, G g){
        generate(m);
    }

    template<class Matrix>
    inline void generate(tiles<Matrix>& m){
        split_d(m);
        int size = m.data.size();
        for(int i = 0; i < size; i++){
            fill_random(m[i]);
        }
    }

    template<class Matrix>
    inline void remove_rows(tiles<Matrix>& m, size_type i, difference_type k){
        assert(false); printf("ERROR: NOT TESTED (BLOCKED REMOVE ROWS)\n");
    }

    template<class Matrix>
    inline void remove_cols(tiles<Matrix>& m, size_type j, difference_type k){
        assert(false); printf("ERROR: NOT TESTED (BLOCKED REMOVE COLS)\n");
    }

    template <class Matrix>
    inline void add_inplace(tiles<Matrix>& lhs, const tiles<Matrix>& rhs){
        split_d(lhs); split_d(rhs);
        int size = lhs.data.size();
        for(int i = 0; i < size; i++){
            lhs[i] += rhs[i];
        }
    }

    template <class Matrix>
    inline void sub_inplace(tiles<Matrix>& lhs, const tiles<Matrix>& rhs){ 
        split_d(lhs); split_d(rhs);
        int size = lhs.data.size();
        for(int i = 0; i < size; i++){
            lhs[i] -= rhs[i];
        }
    }

    template <class MatrixA, class MatrixB>
    inline void mul_inplace(tiles<MatrixA>& m, const tiles<MatrixB>& rhs){
        assert(false); printf("ERROR: NOT TESTED (GEMM BLOCKED INPLACE)\n");
    }

    template <class Matrix>
    inline void mul_inplace(tiles<Matrix>& a, const scalar_type& rhs) { 
        split_d(a);
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            a[i] *= rhs;
        }
    }

    template <class Matrix>
    inline void div_inplace(tiles<Matrix>& a, const scalar_type& rhs){
        split_d(a);
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            a[i] /= rhs;
        }
    }

    template <class Matrix>
    std::ostream& operator << (std::ostream& o, tiles<Matrix> const& m){
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

    template <class Matrix>
    inline const tiles<Matrix> operator * (tiles<Matrix> lhs, const tiles<Matrix>& rhs){ 
        return (lhs *= rhs); 
    }

    template<class Matrix, typename T> 
    inline const tiles<Matrix> operator * (tiles<Matrix> lhs, const T& rhs){ 
        return (lhs *= rhs); 
    }

    template<class Matrix, typename T> 
    inline const tiles<Matrix> operator * (const T& lhs, tiles<Matrix> rhs){ 
        return (rhs *= lhs); 
    }

    template<class Matrix> inline size_type num_rows(const tiles<Matrix>& m){ 
        return m.num_rows(); 
    }

    template<class Matrix> inline size_type num_cols(const tiles<Matrix>& m){
        return m.num_cols(); 
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type 
#endif
