#ifndef AMBIENT_NUMERIC_TILES_HPP
#define AMBIENT_NUMERIC_TILES_HPP

#include "ambient/numeric/matrix/tiles.h"
#include "ambient/numeric/matrix/tiles_algorithms.hpp"

// {{{ tiles< subset_view<Matrix> >
#define size_type typename tiles<subset_view<Matrix> >::size_type

namespace ambient { namespace numeric {

    template<class Matrix>
    inline tiles<subset_view<Matrix> > tiles<subset_view<Matrix> >::subset(size_type i, size_type j, size_type mt, size_type nt) const {
        if(mt == 0 || nt == 0){
            tiles<subset_view<Matrix> > s;
            s.mt = s.nt = s.rows = s.cols = 0;
            return s;
        }

        tiles<subset_view<Matrix> > s;
        s.data.reserve(mt*nt);
        s.mt = mt; s.nt = nt;
        s.rows = (mt-1)*AMBIENT_IB + tile(i+mt-1,0).num_rows(); 
        s.cols = (nt-1)*AMBIENT_IB + tile(0,j+nt-1).num_cols();

        for(int jj = j; jj < j + nt; jj++)
        for(int ii = i; ii < i + mt; ii++)
        s.data.push_back(tile(ii,jj));
        return s;
    }

    template<class Matrix>
    template<class MatrixB>
    inline tiles<subset_view<Matrix> >& tiles<subset_view<Matrix> >::operator += (const tiles<MatrixB>& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix>
    template<class MatrixB>
    inline tiles<subset_view<Matrix> >& tiles<subset_view<Matrix> >::operator -= (const tiles<MatrixB>& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix>
    inline Matrix& tiles<subset_view<Matrix> >::tile(size_type i, size_type j){
        return (Matrix&)this->data[i + mt*j];
    }

    template<class Matrix>
    inline const Matrix& tiles<subset_view<Matrix> >::tile(size_type i, size_type j) const {
        return (Matrix&)this->data[i + mt*j];
    }

    template<class Matrix>
    inline Matrix& tiles<subset_view<Matrix> >::operator[](size_type k){
        return (Matrix&)this->data[k];
    }

    template<class Matrix>
    inline const Matrix& tiles<subset_view<Matrix> >::operator[](size_type k) const {
        return (Matrix&)this->data[k];
    }

    template<class Matrix>
    inline size_type tiles<subset_view<Matrix> >::num_rows() const {
        return this->rows;
    }

    template<class Matrix>
    inline size_type tiles<subset_view<Matrix> >::num_cols() const {
        return this->cols;
    }

    template<class Matrix>
    inline Matrix& tiles<subset_view<Matrix> >::locate(size_type i, size_type j){
        return this->tile(i/AMBIENT_IB, j/AMBIENT_IB);
    }

    template<class Matrix>
    inline const Matrix& tiles<subset_view<Matrix> >::locate(size_type i, size_type j) const {
        return this->tile(i/AMBIENT_IB, j/AMBIENT_IB);
    }

    template<class Matrix>
    inline size_t tiles<subset_view<Matrix> >::addr(size_type i, size_type j) const {
        return locate(i,j).addr(i % AMBIENT_IB, j % AMBIENT_IB);
    }

} }

#undef size_type
// }}}

// {{{ tiles< Matrix >
#define size_type   typename tiles<Matrix>::size_type
#define value_type  typename tiles<Matrix>::value_type
#define scalar_type typename tiles<Matrix>::scalar_type

namespace ambient { namespace numeric {

    template<class Matrix>
    inline void* tiles<Matrix>::operator new (size_t size){
        return ambient::pool.malloc<tiles<Matrix> >();
    }

    template<class Matrix>
    inline void tiles<Matrix>::operator delete (void* ptr){
        ambient::pool.free<tiles<Matrix> >(ptr); 
    }

    template <class Matrix>
    inline tiles<Matrix> tiles<Matrix>::identity_matrix(size_type size){
        tiles t(size, size);
        int tailn = __a_mod(size, AMBIENT_IB);
        for(int i = 0; i < t.nt-1; i++)
            fill_identity(*t.data[i*t.nt + i]);
        fill_identity(*t.data[t.nt*t.nt-1]);
        return t;
    }

    template <class Matrix>
    inline tiles<Matrix>* tiles<Matrix>::new_identity_matrix(size_type size){
        tiles* t = new tiles(size, size);
        int tailn = __a_mod(size, AMBIENT_IB);
        for(int i = 0; i < t->nt-1; i++)
            fill_identity(*t->data[i*t->nt + i]);
        fill_identity(*t->data[t->nt*t->nt-1]);
        return t;
    }

    template <class Matrix>
    inline tiles<Matrix>::tiles()
    : rows(0), cols(0), mt(0), nt(0)
    {
    }

    template <class Matrix>
    inline tiles<Matrix>::tiles(Matrix* a)
    : rows(a->num_rows()), cols(a->num_cols()), mt(__a_ceil(rows/AMBIENT_IB)), nt(__a_ceil(cols/AMBIENT_IB))
    {
        this->data.push_back(a);
    }

    template <class Matrix>
    inline tiles<Matrix>::tiles(size_type rows, size_type cols, value_type init_value)
    : rows(rows), cols(cols), mt(__a_ceil(rows/AMBIENT_IB)), nt(__a_ceil(cols/AMBIENT_IB))
    {
        int tailn = __a_mod(cols, AMBIENT_IB);
        int tailm = __a_mod(rows, AMBIENT_IB);
        this->data.reserve(mt*nt);
        
        for(int j = 1; j < nt; j++){
            for(int i = 1; i < mt; i++) 
                this->data.push_back(new Matrix(AMBIENT_IB, AMBIENT_IB, init_value));
            this->data.push_back(new Matrix(tailm, AMBIENT_IB, init_value));
        }
        for(int i = 1; i < mt; i++) 
            this->data.push_back(new Matrix(AMBIENT_IB, tailn, init_value));
        this->data.push_back(new Matrix(tailm, tailn, init_value));
    }

    template <class Matrix>
    inline tiles<subset_view<Matrix> > tiles<Matrix>::subset(size_type i, size_type j, size_type mt, size_type nt) const {
        if(mt == 0 || nt == 0){
            tiles<subset_view<Matrix> > s;
            s.mt = s.nt = s.rows = s.cols = 0;
            return s;
        }

        tiles<subset_view<Matrix> > s;
        s.data.reserve(mt*nt);
        s.mt = mt; 
        s.nt = nt;
        s.rows = (mt-1)*AMBIENT_IB + tile(i+mt-1,0).num_rows(); 
        s.cols = (nt-1)*AMBIENT_IB + tile(0,j+nt-1).num_cols();
        for(int jj = j; jj < j + nt; jj++)
        for(int ii = i; ii < i + mt; ii++)
        s.data.push_back(subset_view<Matrix>(tile(ii,jj)));

        return s;
    }

    template <class Matrix>
    inline tiles<Matrix>::tiles(const tiles& a)
    : rows(a.rows), cols(a.cols), mt(a.mt), nt(a.nt)
    {
        int nb = a.data.size();
        for(int k = 0; k < nb; k++) this->data.push_back(new Matrix(a[k]));
    }
    
    template <class Matrix>
    tiles<Matrix>& tiles<Matrix>::operator = (const tiles& rhs){
        tiles c(rhs);
        this->swap(c);
        return *this;
    }

#if 0
    template <class Matrix>
    inline tiles<Matrix>::tiles(tiles&& a){
        this->swap(a);
    }

    template <class Matrix>
    tiles<Matrix>& tiles<Matrix>::operator = (tiles&& rhs){
        this->swap(rhs);
        return *this;
    }
#endif

    template <class Matrix>
    tiles<Matrix>::~tiles(){
        int size = this->data.size();
        for(int i = 0; i < size; i++) 
            delete this->data[i];
    }

    template<class Matrix>
    inline size_type tiles<Matrix>::num_rows() const {
        return this->rows;
    }

    template<class Matrix>
    inline size_type tiles<Matrix>::num_cols() const {
        return this->cols;
    }

    template<class Matrix>
    inline scalar_type tiles<Matrix>::trace() const {
        return trace(*this);
    }

    template<class Matrix>
    inline void tiles<Matrix>::transpose(){
        transpose_inplace(*this);
    }

    template<class Matrix>
    inline void tiles<Matrix>::conj(){
        conj_inplace(*this);
    }

    template<class Matrix>
    inline bool tiles<Matrix>::empty() const {
        return (this->rows == 0); 
    }

    template<class Matrix>
    inline void tiles<Matrix>::swap(tiles& r){
        std::swap(r.data,      this->data);
        std::swap(r.rows,      this->rows);
        std::swap(r.cols,      this->cols);
        std::swap(r.mt,        this->mt);
        std::swap(r.nt,        this->nt);
    }

    template<class Matrix>
    inline void tiles<Matrix>::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<class Matrix>
    inline void tiles<Matrix>::remove_rows(size_type i, size_type k){
        remove_rows(*this, i, k);
    }

    template<class Matrix>
    inline void tiles<Matrix>::remove_cols(size_type j, size_type k){
        remove_cols(*this, j, k); 
    }

    template<class Matrix, class OtherMatrix>
    bool operator == (const tiles<Matrix>& a, const tiles<OtherMatrix>& b){
        bool result = true;
        if(a.data.size() != b.data.size()){
            printf("Blocks are different: %lu against %lu\n", a.data.size(), b.data.size());
            return false;
        }
        int size = a.data.size();
        for(int i = 0; i < size; i++){
            if(a[i] == b[i]) continue;
            else result = false;
        }
        return result;
    }

    template<class Matrix>
    inline Matrix& tiles<Matrix>::tile(size_type i, size_type j){
        return (*this)[i + mt*j];
    }

    template<class Matrix>
    inline const Matrix& tiles<Matrix>::tile(size_type i, size_type j) const {
        return (*this)[i + mt*j];
    }

    template<class Matrix>
    inline Matrix& tiles<Matrix>::locate(size_type i, size_type j){
        return this->tile(i/AMBIENT_IB, j/AMBIENT_IB);
    }

    template<class Matrix>
    inline const Matrix& tiles<Matrix>::locate(size_type i, size_type j) const {
        return this->tile(i/AMBIENT_IB, j/AMBIENT_IB);
    }

    template<class Matrix>
    inline size_t tiles<Matrix>::addr(size_type i, size_type j) const {
        return locate(i,j).addr(i % AMBIENT_IB, j % AMBIENT_IB);
    }

    template<class Matrix>
    inline Matrix& tiles<Matrix>::operator[] (size_type k){
        assert(k < data.size());
        return *this->data[k];
    }

    template<class Matrix>
    inline const Matrix& tiles<Matrix>::operator[] (size_type k) const {
        assert(k < data.size());
        return *this->data[k];
    }

    template<class Matrix>
    template<class MatrixB>
    inline tiles<Matrix>::operator tiles<MatrixB> () const {
        tiles<MatrixB> c;
        std::vector<MatrixB*> data;
        data.reserve(this->mt*this->nt);
        for(int j = 0; j < this->nt; j++)
            for(int i = 0; i < this->mt; i++)
                data.push_back(new MatrixB((MatrixB)this->tile(i,j)));
        std::swap(c.data, data);
        c.rows = this->num_rows();
        c.cols = this->num_cols();
        c.mt   = this->mt;
        c.nt   = this->nt;
        return c;
    }

    template<class Matrix>
    template<class MatrixB>
    inline tiles<Matrix>& tiles<Matrix>::operator += (const tiles<MatrixB>& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix>
    template<class MatrixB>
    inline tiles<Matrix>& tiles<Matrix>::operator -= (const tiles<MatrixB>& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix>
    template <typename T2> 
    inline tiles<Matrix>& tiles<Matrix>::operator *= (const T2& t){
        mul_inplace(*this, t);
        return *this;
    }

    template<class Matrix>
    template <typename T2> 
    inline tiles<Matrix>& tiles<Matrix>::operator /= (const T2& t){
        div_inplace(*this, t);
        return *this;
    }

    template<class Matrix>
    inline value_type& tiles<Matrix>::operator() (size_type i, size_type j){
        size_type tiley = this->data[0]->num_rows();
        size_type tilex = this->data[0]->num_cols();
        int mb = __a_ceil(rows/tiley);
        return (*this->data[mb*(int)(j/tilex) + (int)(i/tiley)])(i % tiley, j % tilex);
    }

    template<class Matrix>
    inline const value_type& tiles<Matrix>::operator() (size_type i, size_type j) const {
        size_type tiley = this->data[0]->num_rows();
        size_type tilex = this->data[0]->num_cols();
        int mb = __a_ceil(rows/tiley);
        return (*this->data[mb*(int)(j/tilex) + (int)(i/tiley)])(i % tiley, j % tilex);
    }

} }

#undef size_type
#undef value_type
#undef scalar_type
// }}}

// {{{ tiles< diagonal_matrix<T> >
#define size_type   typename tiles<diagonal_matrix<T> >::size_type
#define value_type  typename tiles<diagonal_matrix<T> >::value_type
#define scalar_type typename tiles<diagonal_matrix<T> >::scalar_type

namespace ambient { namespace numeric {

    template<typename T>
    inline void* tiles<diagonal_matrix<T> >::operator new (size_t size){
        return ambient::pool.malloc<tiles<diagonal_matrix<T> > >();
    }

    template<typename T>
    inline void tiles<diagonal_matrix<T> >::operator delete (void* ptr){
        ambient::pool.free<tiles<diagonal_matrix<T> > >(ptr); 
    }

    template <typename T>
    inline tiles<diagonal_matrix<T> >::tiles(size_type size, value_type init_value)
    : size(size), nt(__a_ceil(size/AMBIENT_IB))
    {
        int tailm = __a_mod(size, AMBIENT_IB);
        this->data.reserve(nt);
        
        for(int i = 1; i < nt; i++) 
            this->data.push_back(new diagonal_matrix<T>(AMBIENT_IB, init_value));
        this->data.push_back(new diagonal_matrix<T>(tailm, init_value));
    }

    template <typename T>
    inline tiles<diagonal_matrix<T> >::tiles(const tiles& a)
    : size(a.size), nt(a.nt)
    {
        int nb = a.data.size();
        for(int k = 0; k < nb; k++) this->data.push_back(new diagonal_matrix<T>(a[k]));
    }
    
    template <typename T>
    tiles<diagonal_matrix<T> >& tiles<diagonal_matrix<T> >::operator = (const tiles& rhs){
        tiles c(rhs);
        this->swap(c);
        return *this;
    }

#if 0
    template <typename T>
    inline tiles<diagonal_matrix<T> >::tiles(tiles&& a){
        this->swap(a);
    }

    template <typename T>
    tiles<diagonal_matrix<T> >& tiles<diagonal_matrix<T> >::operator = (tiles&& rhs){
        this->swap(rhs);
        return *this;
    }
#endif

    template <typename T>
    tiles<diagonal_matrix<T> >::~tiles(){
        int size = this->data.size();
        for(int i = 0; i < size; i++) delete this->data[i];
    }

    template<typename T>
    inline size_type tiles<diagonal_matrix<T> >::num_rows() const {
        return this->size;
    }

    template<typename T>
    inline size_type tiles<diagonal_matrix<T> >::num_cols() const {
        return this->size;
    }

    template<typename T>
    inline void tiles<diagonal_matrix<T> >::swap(tiles& r){
        std::swap(r.data,   this->data);
        std::swap(r.size,   this->size);
        std::swap(r.nt,     this->nt);
    }

    template<typename T>
    inline void tiles<diagonal_matrix<T> >::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<typename T>
    inline void tiles<diagonal_matrix<T> >::remove_rows(size_type i, size_type k){
        remove_rows(*this, i, k);
    }

    template<typename T>
    inline void tiles<diagonal_matrix<T> >::remove_cols(size_type j, size_type k){
        remove_cols(*this, j, k); 
    }

    template<typename T>
    inline diagonal_matrix<T>& tiles<diagonal_matrix<T> >::operator[] (size_type k){
        return *this->data[k];
    }

    template<typename T>
    inline const diagonal_matrix<T>& tiles<diagonal_matrix<T> >::operator[] (size_type k) const {
        return *this->data[k];
    }

    template<typename T>
    inline value_type& tiles<diagonal_matrix<T> >::operator() (size_type i, size_type j){
        size_type tile = this->data[0]->num_rows();
        return (*this->data[(int)(i/tile)])(i % tile, i % tile);
    }

    template<typename T>
    inline const value_type& tiles<diagonal_matrix<T> >::operator() (size_type i, size_type j) const {
        size_type tile = this->data[0]->num_rows();
        return (*this->data[(int)(i/tile)])(i % tile, i % tile);
    }

} }

#undef size_type
#undef value_type
#undef scalar_type
// }}}

#endif
