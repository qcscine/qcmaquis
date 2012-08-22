#ifndef __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__
#define __AMBIENT_NUMERIC_TILES_ALGORITHMS_HPP__

#include "ambient/numeric/matrix/tiles.h"

#define value_type      typename tiles<Matrix>::value_type
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
inline void __a_reduce(std::vector<T*>& seq){
    if(seq.size() == 1) return;
    for(int stride = 1; stride < seq.size(); stride *= 2)
        for(int k = stride; k < seq.size(); k += stride*2)
            *seq[k-stride] += *seq[k];
}

namespace ambient { namespace numeric {

    template<class Matrix>
    bool is_hermitian(const tiles<Matrix>& m){
        return false;
    }

    template<class Matrix>
    inline void merge(tiles<Matrix>& a){
        if(a.data.size() == 1) return;
        int mb = __a_ceil(a.rows/TILE_SIZE);
        int nb = __a_ceil(a.cols/TILE_SIZE);

        std::vector<Matrix*> merge; 
        merge.push_back(new Matrix(a.rows, a.cols));

        for(int j = 0; j < nb; j++){
            for(int i = 0; i < mb; i++){
                Matrix* src = a.data[i+j*mb];
                copy(*merge[0], i*TILE_SIZE, j*TILE_SIZE, *src, 0, 0, src->num_rows(), src->num_cols());
                delete src;
            }
        }
        std::swap(a.data, merge);
    }

    template<class Matrix>
    inline void split(tiles<Matrix>& a){
        if(a.data.size() != 1) return;
        int mb = __a_ceil(a.rows/TILE_SIZE);
        int nb = __a_ceil(a.cols/TILE_SIZE);

        tiles<Matrix> split(a.rows, a.cols);
        Matrix& src = a[0];

        for(int j = 0; j < nb; j++){
            for(int i = 0; i < mb; i++){
                Matrix& dst = split[i+j*mb];
                copy(dst, 0, 0, src, i*TILE_SIZE, j*TILE_SIZE, dst.num_rows(), dst.num_cols());
            }
        }
        swap(a, split);
    }

    template<class MatrixA, class MatrixB, class MatrixC>
    inline void gemm(const tiles<MatrixA>& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        int mb = __a_ceil(c.rows/TILE_SIZE);
        int nb = __a_ceil(c.cols/TILE_SIZE);
        int kb = __a_ceil(a.cols/TILE_SIZE);
        int lda = __a_ceil(a.rows/TILE_SIZE);
        int ldb = __a_ceil(b.rows/TILE_SIZE);

        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                std::vector<MatrixC*> ctree; ctree.reserve(kb);
                ctree.push_back(&c[i + mb*j]);
                size_t rows = c[i + mb*j].num_rows();
                size_t cols = c[i + mb*j].num_cols();
                for(int k = 1; k < kb; k++) 
                    ctree.push_back(new MatrixC(rows, cols));
                for(int k = 0; k < kb; k++){
                    const MatrixA& ab = a[i + k*lda];
                    const MatrixB& bb = b[k + j*ldb];
                    gemm(ab, bb, *ctree[k]);
                }
                __a_reduce(ctree);
                for(int k = 1; k < kb; k++) 
                    delete ctree[k];
            }
        }
    }

    template<class MatrixA, class MatrixC, typename T>
    inline void gemm(const tiles<MatrixA>& a, const tiles<diagonal_matrix<T> >& b, tiles<MatrixC>& c){
        int mb = __a_ceil(c.rows/TILE_SIZE);
        int nb = __a_ceil(c.cols/TILE_SIZE);
        int lda = __a_ceil(a.rows/TILE_SIZE);

        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                gemm(a[i + j*lda], b[j], c[i + mb*j]);
            }
        }
    }

    template<class MatrixB, class MatrixC, typename T>
    inline void gemm(const tiles<diagonal_matrix<T> >& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        int mb = __a_ceil(c.rows/TILE_SIZE);
        int nb = __a_ceil(c.cols/TILE_SIZE);
        int ldb = __a_ceil(b.rows/TILE_SIZE);

        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                gemm(a[i], b[i + j*ldb], c[i + mb*j]);
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
    inline void syev(const tiles<Matrix>& a, tiles<Matrix>& evecs, tiles<DiagonalMatrix>& evals){
        syev(a[0], evecs[0], evals[0]);
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
    inline tiles<Matrix> exp(const tiles<Matrix>& m, const value_type& alfa = 1.){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
    }

    template<class Matrix> 
    inline void resize(tiles<Matrix>& m, size_t rows, size_t cols){ 
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
        swap(m, r);
    }

    template<class Matrix> 
    inline void sqrt_inplace(tiles<Matrix>& m){
        assert(false); printf("ERROR: NOT TESTED (SQRT DIAG)\n");
    }

    template<class Matrix> 
    inline scalar_type trace(const tiles<Matrix>& m){
        int nb = __a_ceil(m.num_cols()/TILE_SIZE);
        int mb = __a_ceil(m.num_rows()/TILE_SIZE);
        int size = std::min(nb, mb);
        scalar_type result(0.);
        std::vector<scalar_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts.push_back(trace(m[k + k*mb]));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline real_type norm_square(const tiles<Matrix>& m){ 
        int size = __a_ceil(m.num_cols()/TILE_SIZE) * __a_ceil(m.num_rows()/TILE_SIZE);
        real_type result(0.);
        std::vector<real_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts.push_back(norm_square(m[k]));
        for(int k = 0; k < size; k++) result += parts[k];
        return result;
    }

    template <class Matrix>
    inline scalar_type overlap(const tiles<Matrix>& a, const tiles<Matrix>& b){
        int size = __a_ceil(a.num_cols()/TILE_SIZE) * __a_ceil(a.num_rows()/TILE_SIZE);
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
        int nb = __a_ceil(a.num_cols()/TILE_SIZE);
        int mb = __a_ceil(a.num_rows()/TILE_SIZE);
        std::vector<Matrix*> t;
        t.reserve(mb*nb);
        for(int i = 0; i < mb; i++){
            for(int j = 0; j < nb; j++){
                Matrix& block = a[i+mb*j];
                transpose_inplace(block);
                t.push_back(&block);
            }
        }
        std::swap(a.rows, a.cols);
        std::swap(a.data, t);
    }

    template<class Matrix>
    inline tiles<transpose_view<Matrix> > transpose(const tiles<Matrix>& a){
        tiles<transpose_view<Matrix> > t;
        int nb = __a_ceil(a.num_cols()/TILE_SIZE);
        int mb = __a_ceil(a.num_rows()/TILE_SIZE);
        std::vector<transpose_view<Matrix>*> data;
        data.reserve(mb*nb);
        for(int i = 0; i < mb; i++)
            for(int j = 0; j < nb; j++)
                data.push_back(new transpose_view<Matrix>(a[i+mb*j]));
        std::swap(t.data, data);
        t.cols = a.num_rows();
        t.rows = a.num_cols();
        return t;
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
