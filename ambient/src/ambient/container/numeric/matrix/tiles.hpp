/*
 * Copyright Institute for Theoretical Physics, ETH Zurich 2015.
 * Distributed under the Boost Software License, Version 1.0.
 *
 * Permission is hereby granted, free of charge, to any person or organization
 * obtaining a copy of the software and accompanying documentation covered by
 * this license (the "Software") to use, reproduce, display, distribute,
 * execute, and transmit the Software, and to prepare derivative works of the
 * Software, and to permit third-parties to whom the Software is furnished to
 * do so, all subject to the following:
 *
 * The copyright notices in the Software and this entire statement, including
 * the above license grant, this restriction and the following disclaimer,
 * must be included in all copies of the Software, in whole or in part, and
 * all derivative works of the Software, unless such copies or derivative
 * works are solely in the form of machine-executable object code generated by
 * a source language processor.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT
 * SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE
 * FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE,
 * ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER
 * DEALINGS IN THE SOFTWARE.
 */

#ifndef AMBIENT_CONTAINER_NUMERIC_TILES_HPP
#define AMBIENT_CONTAINER_NUMERIC_TILES_HPP

#include "ambient/container/numeric/matrix/tiles.h"
#include "ambient/container/numeric/matrix/tiles_algorithms.hpp"

// {{{ tiles< subset_view<Matrix> >
#define size_type typename tiles<subset_view<Matrix>, IB>::size_type

namespace ambient { namespace numeric {

    template<class Matrix, int IB>
    inline tiles<subset_view<Matrix>, IB> tiles<subset_view<Matrix>, IB>::subset(size_type i, size_type j, size_type mt, size_type nt) const {
        if(mt == 0 || nt == 0){
            tiles<subset_view<Matrix>, IB> s;
            s.mt = s.nt = s.rows = s.cols = 0;
            return s;
        }

        tiles<subset_view<Matrix>, IB> s;
        s.data.reserve(mt*nt);
        s.mt = mt; s.nt = nt;
        s.rows = (mt-1)*IB + tile(i+mt-1,0).num_rows(); 
        s.cols = (nt-1)*IB + tile(0,j+nt-1).num_cols();

        for(int jj = j; jj < j + nt; jj++)
        for(int ii = i; ii < i + mt; ii++)
        s.data.push_back(tile(ii,jj));
        return s;
    }

    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<subset_view<Matrix>, IB>& tiles<subset_view<Matrix>, IB>::operator += (const tiles<MatrixB, IB>& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<subset_view<Matrix>, IB>& tiles<subset_view<Matrix>, IB>::operator -= (const tiles<MatrixB, IB>& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }
    
    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<subset_view<Matrix>, IB>& tiles<subset_view<Matrix>, IB>::operator = (const tiles<MatrixB, IB>& rhs){
        int size = rhs.data.size();
        for(int i = 0; i < size; i++)
            (*this)[i] = rhs[i];
        return *this;
    }
    
    template<class Matrix, int IB>
    inline tiles<subset_view<Matrix>, IB>& tiles<subset_view<Matrix>, IB>::operator = (const tiles<subset_view<Matrix>, IB>& rhs){
        int size = rhs.data.size();
        for(int i = 0; i < size; i++)
            (*this)[i] = rhs[i];
        return *this;
    }

    template<class Matrix, int IB>
    inline Matrix& tiles<subset_view<Matrix>, IB>::tile(size_type i, size_type j){
        return (Matrix&)this->data[i + mt*j];
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<subset_view<Matrix>, IB>::tile(size_type i, size_type j) const {
        return (Matrix&)this->data[i + mt*j];
    }

    template<class Matrix, int IB>
    inline Matrix& tiles<subset_view<Matrix>, IB>::operator[](size_type k){
        return (Matrix&)this->data[k];
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<subset_view<Matrix>, IB>::operator[](size_type k) const {
        return (Matrix&)this->data[k];
    }

    template<class Matrix, int IB>
    inline size_type tiles<subset_view<Matrix>, IB>::num_rows() const {
        return this->rows;
    }

    template<class Matrix, int IB>
    inline size_type tiles<subset_view<Matrix>, IB>::num_cols() const {
        return this->cols;
    }

    template<class Matrix, int IB>
    inline Matrix& tiles<subset_view<Matrix>, IB>::locate(size_type i, size_type j){
        return this->tile(i/IB, j/IB);
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<subset_view<Matrix>, IB>::locate(size_type i, size_type j) const {
        return this->tile(i/IB, j/IB);
    }

    template<class Matrix, int IB>
    inline size_t tiles<subset_view<Matrix>, IB>::addr(size_type i, size_type j) const {
        return locate(i,j).addr(i % IB, j % IB);
    }

} }

#undef size_type
// }}}

// {{{ tiles< Matrix >
#define size_type   typename tiles<Matrix,IB>::size_type
#define value_type  typename tiles<Matrix,IB>::value_type
#define scalar_type typename tiles<Matrix,IB>::scalar_type

namespace ambient { namespace numeric {

    template <class Matrix, int IB>
    inline tiles<Matrix, IB> tiles<Matrix, IB>::identity_matrix(size_type size){
        tiles t(size, size);
        for(int i = 0; i < t.nt; i++)
            fill_identity(*t.data[i*t.nt + i]);
        return t;
    }

    template <class Matrix, int IB>
    inline tiles<Matrix, IB>::tiles()
    : rows(0), cols(0), mt(0), nt(0)
    {
    }

    template <class Matrix, int IB>
    inline tiles<Matrix, IB>::tiles(Matrix* a)
    : rows(a->num_rows()), cols(a->num_cols()), mt(__a_ceil(rows/IB)), nt(__a_ceil(cols/IB))
    {
        this->data.push_back(a);
    }

    template <class Matrix, int IB>
    inline tiles<Matrix, IB>::tiles(size_type rows, size_type cols, value_type init_value)
    : rows(rows), cols(cols), mt(__a_ceil(rows/IB)), nt(__a_ceil(cols/IB))
    {
        int tailn = __a_mod(cols, IB);
        int tailm = __a_mod(rows, IB);
        this->data.reserve(mt*nt);
        
        for(int j = 1; j < nt; j++){
            for(int i = 1; i < mt; i++)
                this->data.push_back(new Matrix(IB, IB, init_value));
            this->data.push_back(new Matrix(tailm, IB, init_value));
        }
        for(int i = 1; i < mt; i++) 
            this->data.push_back(new Matrix(IB, tailn, init_value));
        this->data.push_back(new Matrix(tailm, tailn, init_value));
    }

    template <class Matrix, int IB>
    inline tiles<subset_view<Matrix>, IB> tiles<Matrix, IB>::subset(size_type i, size_type j, size_type mt, size_type nt) const {
        if(mt == 0 || nt == 0){
            tiles<subset_view<Matrix>, IB> s;
            s.mt = s.nt = s.rows = s.cols = 0;
            return s;
        }

        tiles<subset_view<Matrix>, IB> s;
        s.data.reserve(mt*nt);
        s.mt = mt; 
        s.nt = nt;
        s.rows = (mt-1)*IB + tile(i+mt-1,0).num_rows(); 
        s.cols = (nt-1)*IB + tile(0,j+nt-1).num_cols();
        for(int jj = j; jj < j + nt; jj++)
        for(int ii = i; ii < i + mt; ii++)
        s.data.push_back(subset_view<Matrix>(tile(ii,jj)));

        return s;
    }

    template <class Matrix, int IB>
    inline tiles<Matrix, IB>::tiles(const tiles& a)
    : rows(a.rows), cols(a.cols), mt(a.mt), nt(a.nt)
    {
        int nb = a.data.size();
        for(int k = 0; k < nb; k++) this->data.push_back(new Matrix(a[k]));
    }
    
    template <class Matrix, int IB>
    tiles<Matrix, IB>& tiles<Matrix, IB>::operator = (const tiles& rhs){
        tiles c(rhs);
        this->swap(c);
        return *this;
    }
    
    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<Matrix, IB>& tiles<Matrix, IB>::operator = (const tiles<MatrixB, IB>& rhs){
        int size = rhs.data.size();
        for(int i = 0; i < size; i++)
            (*this)[i] = rhs[i];
        return *this;
    }

    template <class Matrix, int IB>
    tiles<Matrix, IB>::~tiles(){
        int size = this->data.size();
        for(int i = 0; i < size; i++) 
            delete this->data[i];
    }

    template<class Matrix, int IB>
    inline size_type tiles<Matrix, IB>::num_rows() const {
        return this->rows;
    }

    template<class Matrix, int IB>
    inline size_type tiles<Matrix, IB>::num_cols() const {
        return this->cols;
    }

    template<class Matrix, int IB>
    inline scalar_type tiles<Matrix, IB>::trace() const {
        return trace(*this);
    }

    template<class Matrix, int IB>
    inline void tiles<Matrix, IB>::transpose(){
        transpose_inplace(*this);
    }

    template<class Matrix, int IB>
    inline void tiles<Matrix, IB>::conj(){
        conj_inplace(*this);
    }

    template<class Matrix, int IB>
    inline bool tiles<Matrix, IB>::empty() const {
        return (this->rows == 0); 
    }

    template<class Matrix, int IB>
    inline void tiles<Matrix, IB>::swap(tiles& r){
        std::swap(r.data, this->data);
        std::swap(r.rows, this->rows);
        std::swap(r.cols, this->cols);
        std::swap(r.mt,   this->mt);
        std::swap(r.nt,   this->nt);
    }

    template<class Matrix, int IB>
    inline void tiles<Matrix, IB>::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<class Matrix, class OtherMatrix, int IB>
    bool operator == (const tiles<Matrix, IB>& a, const tiles<OtherMatrix, IB>& b){
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

    template<class Matrix, int IB>
    inline Matrix& tiles<Matrix, IB>::tile(size_type i, size_type j){
        return (*this)[i + mt*j];
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<Matrix, IB>::tile(size_type i, size_type j) const {
        return (*this)[i + mt*j];
    }

    template<class Matrix, int IB>
    inline Matrix& tiles<Matrix, IB>::locate(size_type i, size_type j){
        return this->tile(i/IB, j/IB);
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<Matrix, IB>::locate(size_type i, size_type j) const {
        return this->tile(i/IB, j/IB);
    }

    template<class Matrix, int IB>
    inline size_t tiles<Matrix, IB>::addr(size_type i, size_type j) const {
        return locate(i,j).addr(i % IB, j % IB);
    }

    template<class Matrix, int IB>
    inline Matrix& tiles<Matrix, IB>::operator[] (size_type k){
        assert(k < data.size());
        return *this->data[k];
    }

    template<class Matrix, int IB>
    inline const Matrix& tiles<Matrix, IB>::operator[] (size_type k) const {
        assert(k < data.size());
        return *this->data[k];
    }

    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<Matrix, IB>::operator tiles<MatrixB, IB> () const {
        tiles<MatrixB, IB> c;
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

    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<Matrix, IB>& tiles<Matrix, IB>::operator += (const tiles<MatrixB, IB>& rhs){
        add_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix, int IB>
    template<class MatrixB>
    inline tiles<Matrix, IB>& tiles<Matrix, IB>::operator -= (const tiles<MatrixB, IB>& rhs){
        sub_inplace(*this, rhs);
        return *this;
    }

    template<class Matrix, int IB>
    template <typename T2> 
    inline tiles<Matrix, IB>& tiles<Matrix, IB>::operator *= (const T2& t){
        mul_inplace(*this, t);
        return *this;
    }

    template<class Matrix, int IB>
    template <typename T2> 
    inline tiles<Matrix, IB>& tiles<Matrix, IB>::operator /= (const T2& t){
        div_inplace(*this, t);
        return *this;
    }

    template<class Matrix, int IB>
    inline value_type& tiles<Matrix, IB>::operator() (size_type i, size_type j){
        size_type tiley = this->data[0]->num_rows();
        size_type tilex = this->data[0]->num_cols();
        int mb = __a_ceil(rows/tiley);
        return ambient::load(*this->data[mb*(int)(j/tilex) + (int)(i/tiley)])(i % tiley, j % tilex);
    }

    template<class Matrix, int IB>
    inline const value_type& tiles<Matrix, IB>::operator() (size_type i, size_type j) const {
        size_type tiley = this->data[0]->num_rows();
        size_type tilex = this->data[0]->num_cols();
        int mb = __a_ceil(rows/tiley);
        return ambient::load(*this->data[mb*(int)(j/tilex) + (int)(i/tiley)])(i % tiley, j % tilex);
    }

} }

#undef size_type
#undef value_type
#undef scalar_type
// }}}

// {{{ tiles< diagonal_matrix<T>, IB >
#define size_type   typename tiles<diagonal_matrix<T>, IB>::size_type
#define value_type  typename tiles<diagonal_matrix<T>, IB>::value_type
#define scalar_type typename tiles<diagonal_matrix<T>, IB>::scalar_type

namespace ambient { namespace numeric {

    template <typename T, int IB>
    inline tiles<diagonal_matrix<T>, IB> tiles<diagonal_matrix<T>, IB>::identity_matrix(size_type size){
        return tiles(size, size, 1.);
    }

    template <typename T, int IB>
    inline tiles<diagonal_matrix<T>, IB>::tiles()
    : size(0), nt(0)
    {
    }

    template <typename T, int IB>
    inline tiles<diagonal_matrix<T>, IB>::tiles(size_type rows, size_type cols, value_type init_value)
    : size(rows), nt(__a_ceil(rows/IB))
    {
        assert(rows == cols);
        int tailm = __a_mod(size, IB);
        this->data.reserve(nt);
        
        for(int i = 1; i < nt; i++) 
            this->data.push_back(new diagonal_matrix<T>(IB, IB, init_value));
        this->data.push_back(new diagonal_matrix<T>(tailm, tailm, init_value));
    }

    template <typename T, int IB>
    inline tiles<diagonal_matrix<T>, IB>::tiles(const tiles& a)
    : size(a.size), nt(a.nt)
    {
        int nb = a.data.size();
        for(int k = 0; k < nb; k++) this->data.push_back(new diagonal_matrix<T>(a[k]));
    }
    
    template <typename T, int IB>
    tiles<diagonal_matrix<T>, IB>& tiles<diagonal_matrix<T>, IB>::operator = (const tiles& rhs){
        tiles c(rhs);
        this->swap(c);
        return *this;
    }

    template <typename T, int IB>
    tiles<diagonal_matrix<T>, IB>::~tiles(){
        int size = this->data.size();
        for(int i = 0; i < size; i++) delete this->data[i];
    }

    template<typename T, int IB>
    inline std::pair<const value_type*,const value_type*> tiles<diagonal_matrix<T>, IB>::diagonal() const {
        return std::make_pair(this->begin(), this->end());
    }
        
    template<typename T, int IB>
    inline const value_type* tiles<diagonal_matrix<T>, IB>::begin() const {
        return &(*this)[0][0];
    }
        
    template<typename T, int IB>
    inline const value_type* tiles<diagonal_matrix<T>, IB>::end() const {
        return (&(*this)[0][0] + (*this)[0].num_rows());
    }

    template<typename T, int IB>
    inline size_type tiles<diagonal_matrix<T>, IB>::num_rows() const {
        return this->size;
    }

    template<typename T, int IB>
    inline size_type tiles<diagonal_matrix<T>, IB>::num_cols() const {
        return this->size;
    }

    template<typename T, int IB>
    inline void tiles<diagonal_matrix<T>, IB>::swap(tiles& r){
        std::swap(r.data,   this->data);
        std::swap(r.size,   this->size);
        std::swap(r.nt,     this->nt);
    }

    template<typename T, int IB>
    inline void tiles<diagonal_matrix<T>, IB>::resize(size_type m, size_type n){
        ambient::numeric::resize(*this, m, n);
    }

    template<typename T, int IB>
    inline diagonal_matrix<T>& tiles<diagonal_matrix<T>, IB>::operator[] (size_type k){
        return *this->data[k];
    }

    template<typename T, int IB>
    inline const diagonal_matrix<T>& tiles<diagonal_matrix<T>, IB>::operator[] (size_type k) const {
        return *this->data[k];
    }

    template<typename T, int IB>
    inline value_type& tiles<diagonal_matrix<T>, IB>::operator() (size_type i, size_type j){
        size_type tile = this->data[0]->num_rows();
        return (*this->data[(int)(i/tile)])(i % tile, i % tile);
    }

    template<typename T, int IB>
    inline const value_type& tiles<diagonal_matrix<T>, IB>::operator() (size_type i, size_type j) const {
        size_type tile = this->data[0]->num_rows();
        return (*this->data[(int)(i/tile)])(i % tile, i % tile);
    }

} }

#undef size_type
#undef value_type
#undef scalar_type
// }}}

#endif
