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
        for(int k = stride; k < seq.size(); k += stride*2)
            *seq[k-stride] += *seq[k];
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
        step = std::min(std::min((SPLIT_SIZE*__a_ceil((first+1)/SPLIT_SIZE) - first), 
                                 (SPLIT_SIZE*__a_ceil((second+1)/SPLIT_SIZE) - second)),
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
        int tailm = __a_mod(a.rows, SPLIT_SIZE);
        int tailn = __a_mod(a.cols, SPLIT_SIZE);
        split.reserve(a.mt*a.nt);
        for(int j = 1; j < a.nt; j++){
            for(int i = 1; i < a.mt; i++) 
                split.push_back(new Matrix(SPLIT_SIZE, SPLIT_SIZE));
            split.push_back(new Matrix(tailm, SPLIT_SIZE));
        }
        for(int i = 1; i < a.mt; i++) 
            split.push_back(new Matrix(SPLIT_SIZE, tailn));
        split.push_back(new Matrix(tailm, tailn));

        const Matrix* src = a.data[0];

        if(!a[0].impl->weak())
        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                Matrix& dst = *split[i+j*a.mt];
                copy(dst, 0, 0, *src, i*SPLIT_SIZE, j*SPLIT_SIZE, dst.num_rows(), dst.num_cols());
            }
        }

        delete src;
        std::swap(const_cast<std::vector<Matrix*>&>(a.data), split);
    }

    template<class Matrix>
    inline void split_d(const tiles<transpose_view<Matrix> >& a){
        if(a.data.size() == 1 && (__a_ceil(a.rows/SPLIT_SIZE) != 1 || __a_ceil(a.cols/SPLIT_SIZE) != 1)) printf("We have to split but we don't :(\n");
    }

    template<typename T>
    inline void split_d(const tiles<diagonal_matrix<T> >& a){
        if(a.data.size() != 1) return;
        if(a.nt == 1) return;

        std::vector<diagonal_matrix<T> *> split;
        int tailm = __a_mod(a.size, SPLIT_SIZE);
        split.reserve(a.nt);
        for(int i = 1; i < a.nt; i++) 
            split.push_back(new diagonal_matrix<T>(SPLIT_SIZE));
        split.push_back(new diagonal_matrix<T>(tailm));

        const diagonal_matrix<T>* src = a.data[0];

        if(!a[0].get_data().impl->weak())
        for(int i = 0; i < a.nt; i++){
            diagonal_matrix<T>& dst = *split[i];
            copy(dst.get_data(), 0, 0, src->get_data(), i*SPLIT_SIZE, 0, dst.num_rows(), 1);
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
                const Matrix* src = a.data[i+j*a.mt];
                if(!src->impl->weak())
                copy(*merge[0], i*SPLIT_SIZE, j*SPLIT_SIZE, *src, 0, 0, src->num_rows(), src->num_cols());
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
        int tailm = __a_mod(a.rows, SPLIT_SIZE);
        int tailn = __a_mod(a.cols, SPLIT_SIZE);
        split.reserve(a.mt*a.nt);
        for(int j = 1; j < a.nt; j++){
            for(int i = 1; i < a.mt; i++) 
                split.push_back(new Matrix(SPLIT_SIZE, SPLIT_SIZE));
            split.push_back(new Matrix(tailm, SPLIT_SIZE));
        }
        for(int i = 1; i < a.mt; i++) 
            split.push_back(new Matrix(SPLIT_SIZE, tailn));
        split.push_back(new Matrix(tailm, tailn));

        const Matrix* src = a.data[0];

        if(!a[0].impl->weak())
        for(int j = 0; j < a.nt; j++){
            for(int i = 0; i < a.mt; i++){
                Matrix& dst = *split[i+j*a.mt];
                copy(dst, 0, 0, *src, i*SPLIT_SIZE, j*SPLIT_SIZE, dst.num_rows(), dst.num_cols());
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
            copy(merge[0]->get_data(), i*SPLIT_SIZE, 0, src->get_data(), 0, 0, src->num_rows(), 1);
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
        int tailm = __a_mod(a.size, SPLIT_SIZE);
        split.reserve(a.nt);
        for(int i = 1; i < a.nt; i++) 
            split.push_back(new diagonal_matrix<T>(SPLIT_SIZE));
        split.push_back(new diagonal_matrix<T>(tailm));

        const diagonal_matrix<T>* src = a.data[0];

        if(!a[0].get_data().impl->weak())
        for(int i = 0; i < a.nt; i++){
            diagonal_matrix<T>& dst = *split[i];
            copy(dst.get_data(), 0, 0, src->get_data(), i*SPLIT_SIZE, 0, dst.num_rows(), 1);
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
        copy(out.locate(row.first, col.first), row.first%SPLIT_SIZE, col.first%SPLIT_SIZE, 
             in.locate(row.second, col.second), row.second%SPLIT_SIZE, col.second%SPLIT_SIZE,
             row.step, col.step);
    }

    template<class Matrix>
    inline void copy_s(tiles<Matrix>& out, size_t oi, size_t oj, 
                       const tiles<Matrix>& in, size_t ii, size_t ij,
                       const tiles<Matrix>& alfa, size_t ai, size_t aj,
                       size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= SPLIT_SIZE; aj %= SPLIT_SIZE;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_s(out.locate(row.first, col.first), row.first%SPLIT_SIZE, col.first%SPLIT_SIZE, 
               in.locate(row.second, col.second), row.second%SPLIT_SIZE, col.second%SPLIT_SIZE,
               factor, ai, aj, row.step, col.step);
    }

    template<class Matrix>
    inline void copy_sa(tiles<Matrix>& out, size_t oi, size_t oj, 
                        const tiles<Matrix>& in, size_t ii, size_t ij,
                        const tiles<Matrix>& alfa, size_t ai, size_t aj,
                        size_t m, size_t n)
    {
        const Matrix& factor = alfa.locate(ai, aj); 
        ai %= SPLIT_SIZE; aj %= SPLIT_SIZE;
        for(cross_iterator row(oi,ii,m); !row.end(); ++row)
        for(cross_iterator col(oj,ij,n); !col.end(); ++col)
        copy_sa(out.locate(row.first, col.first), row.first%SPLIT_SIZE, col.first%SPLIT_SIZE, 
                in.locate(row.second, col.second), row.second%SPLIT_SIZE, col.second%SPLIT_SIZE,
                factor, ai, aj, row.step, col.step);
    }

    template<class MatrixA, class MatrixB, class MatrixC>
    inline void gemm(const tiles<MatrixA>& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        split_d(a); split_d(b); split_d(c);

        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                std::vector<MatrixC*> ctree; ctree.reserve(a.nt);
                ctree.push_back(&c[i + c.mt*j]);
                size_t rows = c[i + c.mt*j].num_rows();
                size_t cols = c[i + c.mt*j].num_cols();
                for(int k = 1; k < a.nt; k++) 
                    ctree.push_back(new MatrixC(rows, cols));
                for(int k = 0; k < a.nt; k++){
                    const MatrixA& ab = a[i + k*a.mt];
                    const MatrixB& bb = b[k + j*b.mt];
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
                gemm(a[i + j*a.mt], b[j], c[i + c.mt*j]);
            }
        }
    }

    template<class MatrixB, class MatrixC, typename T>
    inline void gemm(const tiles<diagonal_matrix<T> >& a, const tiles<MatrixB>& b, tiles<MatrixC>& c){
        split_d(a); split_d(b); split_d(c);

        for(int i = 0; i < c.mt; i++){
            for(int j = 0; j < c.nt; j++){
                gemm(a[i], b[i + j*b.mt], c[i + c.mt*j]);
            }
        }
    }

    template<class Matrix, class DiagonalMatrix>
    inline void svd(const tiles<Matrix>& a, tiles<Matrix>& u, tiles<Matrix>& vt, tiles<DiagonalMatrix>& s){
        merge(a); merge(u); merge(vt); merge(s);
        svd(a[0], u[0], vt[0], s[0]);
        split(a); split(u); split(vt); split(s);
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
    inline void qr(const tiles<Matrix>& a, tiles<Matrix>& q, tiles<Matrix>& r){
        merge(a); merge(q); merge(r);
        qr(a[0], q[0], r[0]); 
        split(a); split(q); split(r);
    }

    template<class Matrix> 
    inline void lq(const tiles<Matrix>& a, tiles<Matrix>& l, tiles<Matrix>& q){ 
        merge(a); merge(l); merge(q);
        qr(a[0], l[0], q[0]); 
        split(a); split(l); split(q);
    } 

    template<class Matrix>
    inline tiles<Matrix> exp(const tiles<Matrix>& m, const value_type& alfa = 1.){
        assert(false); printf("ERROR: NOT TESTED (EXP)\n");
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
                r[i+j*r.mt] = m[i+j*m.mt];
            }
        }
        size_t margin = SPLIT_SIZE;
        if(r.mt <= m.mt) margin = __a_mod(rows, SPLIT_SIZE);
        for(int j = 0; j < nb_min; j++){
            Matrix& block = r[mb_min-1 + j*r.mt];
            resize(block, margin, block.num_cols());
        }
        margin = SPLIT_SIZE;
        if(r.nt <= m.nt) margin = __a_mod(cols, SPLIT_SIZE);
        for(int i = 0; i < mb_min; i++){
            Matrix& block = r[i + (nb_min-1)*r.mt];
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
            resize(r[nb_min-1], SPLIT_SIZE, SPLIT_SIZE);
        }else{
            size_t margin = __a_mod(rows, SPLIT_SIZE);
            resize(r[r.nt-1], margin, margin);
        }
        swap(m, r);
    }

    template<class Matrix> 
    inline void sqrt_inplace(tiles<Matrix>& m){
        assert(false); printf("ERROR: NOT TESTED (SQRT DIAG)\n");
    }

    template<class Matrix> 
    inline scalar_type trace(const tiles<Matrix>& m){
        split_d(m);
        int size = std::min(m.nt, m.mt);
        scalar_type result(0.);
        std::vector<scalar_type> parts;
        parts.reserve(size);

        for(int k = 0; k < size; k++) parts.push_back(trace(m[k + k*m.mt]));
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
                Matrix& block = a[i+a.mt*j];
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
        split_d(a);
        tiles<transpose_view<Matrix> > t;
        std::vector<transpose_view<Matrix>*> data;
        data.reserve(a.mt*a.nt);
        for(int i = 0; i < a.mt; i++)
            for(int j = 0; j < a.nt; j++)
                data.push_back(new transpose_view<Matrix>(a[i+a.mt*j]));
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
        split_d(m);
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
