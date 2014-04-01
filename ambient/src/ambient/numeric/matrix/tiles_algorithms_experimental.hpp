/*
 * Ambient, License - Version 1.0 - May 3rd, 2012
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

#ifndef AMBIENT_NUMERIC_TILES_ALGORITHMS_EXPERIMENTAL
#define AMBIENT_NUMERIC_TILES_ALGORITHMS_EXPERIMENTAL

#define value_type      typename tiles<Matrix>::value_type
#define size_type       typename tiles<Matrix>::size_type
#define real_type       typename tiles<Matrix>::real_type
#define scalar_type     typename tiles<Matrix>::scalar_type
#define difference_type typename tiles<Matrix>::difference_type

#define PI_VALUE 3.14159265359 

namespace ambient { namespace numeric {

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
            /* {{{ explicit
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

            split(x);
            split(y);

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
        
        tiles<DiagonalMatrix> tp(k, k);
        tiles<DiagonalMatrix> tq(k, k);
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
        tiles<DiagonalMatrix> e(k,k);

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
        split(s); split(u2); split(v2);

        ambient::numeric::gemm(u1, u2, u);
        ambient::numeric::gemm(v2, v1, v);
    }

    template<int alfa, int beta, class MatrixA, class MatrixB, class MatrixC>
    inline void gemv(const tiles<MatrixA>&  a, int ai, int aj, 
                     const tiles<MatrixB>&  b, int bi, int bj, 
                           MatrixC&  c, int ci, int cj,
                           int m, int n)
    {
        std::cout << " JE SUIS LA ALEX SVD " << std::endl; // Tim
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

    template<int alfa, int beta, class Matrix>
    inline void gemv(const tiles<Matrix>&  a, int ai, int aj,
                     const tiles<Matrix>&  b, int bi, int bj,
                           tiles<Matrix>&  c, int ci, int cj,
                           int m, int n)
    {
        std::cout << " JE SUIS LA TIM SVD " << std::endl; // for Tim SVD only
        if(m == 0 || n == 0) return;

        for(cross_iterator row(ai,ci,m); !row.end(); ++row){
            std::vector<Matrix*> ctree;
            for(cross_iterator col(aj,bi,n); !col.end(); ++col){
                Matrix* part = new Matrix(row.step, 1);
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

    template <class Matrix>
    inline bool test_norm(const tiles<Matrix>& a, double epsilon){
        int num_cols = a.num_cols();
        tiles<matrix<double> > norm(1,num_cols); 
           
        for(int i=0; i<a.nt; ++i)
            for(int j=0; j<a.mt; ++j)
                norm_vector(a.tile(i,j),norm.tile(0,j));                    
         
        int size = norm.data.size();
        for(int i = 0; i < size; ++i)
            sqrt_inplace(norm[i]); 

        std::vector<double> parts;
        parts.reserve(num_cols);

        for(int i = 0; i < size; ++i)
            parts.push_back(max_vector(norm[i]));
        
        std::vector<double>::iterator maxit;

        maxit = std::max_element(parts.begin(), parts.end());
        
        if(*maxit < epsilon/(10*sqrt(2/PI_VALUE))){
            std::cout << " JE STOP " << std::endl;
            return false; // we stop
        }
 
        std::cout << " JE CONTINUE (max,epsilon) " << *maxit << " " << epsilon  << std::endl;
        return true;  // we continue
    }
 
    template<class Matrix, class DiagonalMatrix>
    inline void svd_lowrank(const tiles<Matrix>& a, tiles<Matrix>& u, tiles<Matrix>& v, tiles<DiagonalMatrix>& s, real_type epsilon){
        int r = 10; // to template later maybe
        tiles<Matrix> omega(a.num_rows(),r);
        generate_gaussian(omega);
        omega *= a; // omega -> y, line 2
        int j(0);
        tiles<Matrix> Q(1,1);
        tiles<Matrix> yn(a.num_rows(),1);
        while(test_norm(omega,epsilon) == true){ // line 5
            ++j; // line 6 ^_^
            std::cout << " JE SUIS 0 " << std::endl;
            tiles<Matrix> Id = tiles<Matrix>::identity_matrix(Q.num_rows());
            std::cout << " JE SUIS 1 " << std::endl;
            Id -= Q * adjoint(Q) ;
            std::cout << " JE SUIS 2 " << std::endl;
            std::cout << omega << std::endl;
            gemv<1,1>(Id,0,0, omega,0,j, yn,0,0, 1, a.num_rows()); // y := alpha*A*x + beta*y,
            std::cout << omega << std::endl;
            assert(false);
            std::cout << " JE SUIS 3 " << std::endl;
        }
    }

    template<class Matrix>
    inline void generate_gaussian(tiles<Matrix>& a){
        int size = a.data.size();
        for(int i = 0; i < size; i++)
            fill_gaussian(a[i]);
    }

} }

#undef size_type
#undef real_type
#undef scalar_type
#undef difference_type 
#endif
