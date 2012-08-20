#ifndef __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__
#define __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__

#include "ambient/numeric/matrix/tiles.h"

#define size_type       typename tiles<Matrix>::size_type
#define real_type       typename tiles<Matrix>::real_type
#define scalar_type     typename tiles<Matrix>::scalar_type
#define difference_type typename tiles<Matrix>::difference_type

inline size_t __a_mod(size_t size, size_t tile){
    size_t mask[2] = {(size & (tile-1)), tile}; 
    return mask[!mask[0]];
}

inline size_t __a_mod_classic(size_t size, size_t tile){
    size_t m = size % tile;
    if(m == 0) m = tile;
    return m;
}

template<typename T>
inline T& __a_reduce(std::vector<T>& seq){
    for(int stride = 1; stride < seq.size(); stride *= 2)
        for(int k = stride; k < seq.size(); k += stride*2)
            seq[k-stride] += seq[k];
    return seq[0];
}

namespace ambient { namespace numeric {

    template<class MatrixA, class MatrixB, typename MatrixC>
    inline void gemm(const tiles<MatrixA>& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        int nb = __a_ceil(c.cols/TILE_SIZE);
        int mb = __a_ceil(c.rows/TILE_SIZE);
        int kb = __a_ceil(a.cols/TILE_SIZE);
        int lda = __a_ceil(a.rows/TILE_SIZE);
        int ldb = __a_ceil(b.rows/TILE_SIZE);

        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                std::vector<MatrixC> ctree;
                ctree.reserve(kb);
                for(int k = 0; k < kb; k++){
                    MatrixA ab = a[i + kb*lda];
                    MatrixB bb = b[kb + j*ldb];
                    ctree.push_back(MatrixC(ab.num_rows(), bb.num_cols()));
                    gemm(ab, bb, ctree[k]);
                }
                c[i + mb*j] = __a_reduce(ctree);
            }
        }
    }

    template<class Matrix, class DiagonalMatrix>
    inline void svd(const tiles<Matrix>& a, tiles<Matrix>& u, tiles<Matrix>& vt, tiles<DiagonalMatrix>& s){
        svd(a[0], u[0], vt[0], s[0]);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void heev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        heev(a[0], evecs[0], evals[0]);
    }

    template<class Matrix, class DiagonalMatrix>
    inline void heev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        heev(a, evecs, evals);
    }

    template<class Matrix> 
    inline void qr(const tiles<Matrix>& a, tiles<Matrix>& q, tiles<Matrix>& r){
       qr(a[0], q[0], r[0]); 
    }

    template<class Matrix> 
    inline void lq(const tiles<Matrix>& a, tiles<Matrix>& l, tiles<Matrix>& q){ 
       qr(a[0], l[0], q[0]); 
    } 

    template<class Matrix>
    inline tiles<Matrix> exp(const tiles<Matrix>& m, const T& alfa = 1.){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
    }

    template<class Matrix> inline void resize(tiles<Matrix>& m, size_t rows, size_t cols){ 
        if(m.num_rows() == rows && m.num_cols() == cols) return;
        tiles<Matrix> r(rows, cols);
        int nb = __a_ceil(cols/TILE_SIZE);
        int mb = __a_ceil(rows/TILE_SIZE);
        int mbo = __a_ceil(m.num_rows()/TILE_SIZE);
        int nbo = __a_ceil(m.num_cols()/TILE_SIZE);
        int mb_min = std::min(mb, mbo);
        int nb_min = std::min(nb, nbo);

        for(int j = 0; j < nb_min; j++){
            for(int i = 0; i < mb_min; i++){
                r[i+j*mb] = m[i+j*mbo];
            }
        }
        if(mb > mbo){
            for(int j = 0; j < nb_min; j++){
                Matrix& block = r[mb_min-1 + j*mb];
                resize(block, TILE_SIZE, block.num_cols());
            }
        }else{
            size_t margin = __a_mod(rows, TILE_SIZE);
            for(int j = 0; j < nb; j++){
                Matrix& block = r[mb-1 + j*mb];
                resize(block, margin, block.num_cols());
            }
        }
        if(nb > nbo){
            for(int i = 0; i < mb_min; i++){
                Matrix& block = r[i + (nb_min-1)*mb];
                resize(block, block.num_rows(), TILE_SIZE);
            }
        }else{
            size_t margin = __a_mod(cols, TILE_SIZE);
            for(int i = 0; i < mb; i++){
                Matrix& block = r[i + (nb-1)*mb];
                resize(block, block.num_rows(), margin);
            }
        }
    }

    template<class Matrix> inline scalar_type trace(const tiles<Matrix>& m){
        int nb = __a_ceil(a.num_cols()/TILE_SIZE);
        int mb = __a_ceil(a.num_rows()/TILE_SIZE);
        int size = std::min(nb, mb);
        Matrix::scalar_type result(0);
        std::vector<Matrix::scalar_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts[k + k*mb] = trace(a[k]);
        for(int k = 0; k < size; k++) result += parts[k];
    }

    template <class Matrix>
    inline real_type norm_square(const tiles<Matrix>& a){ 
        int size = __a_ceil(a.num_cols()/TILE_SIZE) * __a_ceil(a.num_rows()/TILE_SIZE);
        Matrix::scalar_type result(0);
        std::vector<Matrix::scalar_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts[k] = norm_square(a[k]);
        for(int k = 0; k < size; k++) result += parts[k];
    }

    template <class Matrix>
    inline scalar_type overlap(const tiles<Matrix>& a, const tiles<Matrix>& b){
        int size = __a_ceil(a.num_cols()/TILE_SIZE) * __a_ceil(a.num_rows()/TILE_SIZE);
        Matrix::scalar_type result(0);
        std::vector<Matrix::scalar_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts[k] = overlap(a[k], b[k]);
        for(int k = 0; k < size; k++) result += parts[k];
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
        int nb = __a_ceil(a.num_cols()/TILE_SIZE);
        int mb = __a_ceil(a.num_rows()/TILE_SIZE);
        std::vector<Matrix*> t;
        t.reserve(mb*nb);
        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                Matrix* block = a[i+mb*j];
                transpose_inplace(block);
                t.push_back(block);
            }
        }
        std::swap(a.rows, a.cols);
        std::swap(a.data, t);
    }

    template<class Matrix>
    inline tiles<Matrix> transpose(const tiles<Matrix>& a){
        tiles<Matrix> t(a.num_cols(), a.num_rows());
        int nb = __a_ceil(a.num_cols()/TILE_SIZE);
        int mb = __a_ceil(a.num_rows()/TILE_SIZE);
        for(int j = 0; j < nb; i++){
            for(int i = 0; i < mb; j++){
                t[j+i*nb] = transpose(a[i+j*mb]);
            }
        }
    }

    template<class Matrix>
    inline void adjoint_inplace(tiles<Matrix>& m){
        transpose_inplace(m);
    }

    template<class Matrix>
    inline const tiles<Matrix>& adjoint(const tiles<Matrix>& m){
        return transpose(m);
    }

    template<class Matrix, class G>
    inline void generate(tiles<Matrix>& m, G g){
        int size = m.data.size();
        for(int i = 0; i < size; i++){
            generate(m[i], g);
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
        int size = lhs.data.size();
        for(int i = 0; i < size; i++){
            lhs[i] += rhs[i];
        }
    }

    template <class Matrix>
    inline void sub_inplace(tiles<Matrix>& lhs, const tiles<Matrix>& rhs){ 
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
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            a[i] *= rhs;
        }
    }

    template <class Matrix>
    inline void div_inplace(tiles<Matrix>& a, const scalar_type& rhs){
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
