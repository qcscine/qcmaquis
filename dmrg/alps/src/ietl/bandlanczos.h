/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2001-2010 by Rene Villiger <rvilliger@smile.ch>,
*                            Prakash Dayal <prakash@comp-phys.org>,
*                            Matthias Troyer <troyer@comp-phys.org>
*
* This software is part of the ALPS libraries, published under the ALPS
* Library License; you can use, redistribute it and/or modify it under
* the terms of the license, either version 1 or (at your option) any later
* version.
* 
* You should have received a copy of the ALPS Library License along with
* the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
* available from http://alps.comp-phys.org/.
*
* THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR 
* IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY, 
* FITNESS FOR A PARTICULAR PURPOSE, TITLE AND NON-INFRINGEMENT. IN NO EVENT 
* SHALL THE COPYRIGHT HOLDERS OR ANYONE DISTRIBUTING THE SOFTWARE BE LIABLE 
* FOR ANY DAMAGES OR OTHER LIABILITY, WHETHER IN CONTRACT, TORT OR OTHERWISE, 
* ARISING FROM, OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER 
* DEALINGS IN THE SOFTWARE.
*
*****************************************************************************/

/* $Id: bandlanczos.h,v 1.14 2003/09/05 08:12:38 troyer Exp $ */

#ifndef IETL_BANDLANCZOS_H
#define IETL_BANDLANCZOS_H

#include <vector>
#include <algorithm>
#include <ietl/fmatrix.h>
#include <ietl/traits.h>
#include <ietl/iteration.h>
#include <boost/numeric/bindings/lapack/driver/syev.hpp>
#include <boost/numeric/bindings/lapack/driver/heev.hpp>

namespace ietl {  
  
  class indexer { // class indexer starts
  public:
    indexer(int pc);
    ~indexer();
    int cnv(int j);
    void next();
    void deflate(int old);
  private:
    std::vector<int> location_;
    std::vector<bool> deflated_;
    int pc_;
    int lastloc_;
    int lastvec_;
  };
  
  template <class MATRIX, class VS>  // class bandlanczos starts
    class bandlanczos {
    public:
    typedef typename VS::vector_type vector_type;
    typedef typename VS::scalar_type scalar_type;
    typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
    typedef std::vector<scalar_type> eigenvec_type;
    
    bandlanczos(const MATRIX& matrix,const VS& vec, const int& p);
    ~bandlanczos();
    
    template <class GEN, class ITER> void calculate_eigenvalues(ITER& iter,
                                                                const GEN& gen);
    template <class GEN, class ITER> void calculate_eigenvectors(ITER& iter,
                                                                 const GEN& gen);
     
    const std::vector<vector_type>& eigenvectors();
    const std::vector<int>& multiplicities();
    const std::vector<magnitude_type>& eigenvalues();
    
    private:
    bool convergence_test(const int dim, int pc, int evs, magnitude_type dep_tol, magnitude_type ghost_tol, bool ghost_discarding, bool low);
    int calc_approx_eigenvalues(const int& dim, int evs, magnitude_type ghost_tol, bool ghost_discarding, bool low); // store evs eigenvalues, if evs==-1, store all of them
    void getev_(const fortran_int_t dim, double w[], fortran_int_t info, double type);
    void getev_(const fortran_int_t dim, float w[], fortran_int_t info, float type);
    void getev_(const fortran_int_t dim, double w[], fortran_int_t info, std::complex<double> type);
    void getev_(const fortran_int_t dim, float w[], fortran_int_t info, std::complex<float> type);
    
    const MATRIX& matrix_;
    VS vecspace_;
    int p_;
    int n_;
    FortranMatrix<scalar_type> tmat, Tjpr;
    mutable std::vector<magnitude_type> eigenvals;
    mutable std::vector<vector_type> reigenvecs;
    mutable std::vector<eigenvec_type> eigenvecs;
    mutable std::vector<int> multiplicity;
    std::vector<vector_type> v;
   };      
  
  // C L A S S :   I N D E X E R /////////////////////////////////////////////  
  indexer::indexer(int pc) : location_(2*pc+1), deflated_(2*pc+1) {
    pc_=pc;
    lastloc_=pc;
    lastvec_=pc;
    for (int i=0;i<(2*pc+1);i++) {  
      location_[i] = i;
      deflated_[i] = false; 
    }        
  }   
  indexer::~indexer() {}
  
  int indexer::cnv(int j) {      
    return std::find(location_.begin(),location_.end(),j) - location_.begin();
  }
  
  void indexer::next(){
    lastvec_++;
    do {
      lastloc_++;
      if (lastloc_ == location_.size() )
        lastloc_=0;
    } while (deflated_[lastloc_]);
    location_[lastloc_] = lastvec_;
  }
  
  void indexer::deflate(int old) {
    pc_--;
    int place = std::find(location_.begin(),location_.end(),old) - location_.begin();
    deflated_[place] = true;
  }  
  // E N D   O F   C L A S S   I N D E X E R /////////////////////////////////   
     
   // C L A S S :   B A N D L A N C Z O S /////////////////////////////////////   
  template <class MATRIX, class VS>
    bandlanczos<MATRIX, VS>::bandlanczos(const MATRIX& matrix,const VS& vec,
                                         const int& p) : matrix_(matrix), p_(p), 
    vecspace_(vec),v(2*p+1),tmat(1,1),Tjpr(1,1) {
    n_ = vecspace_.vec_dimension();
  }      
  template <class MATRIX, class VS>
    bandlanczos<MATRIX, VS>::~bandlanczos() { }
  
  template <class MATRIX, class VS>
    const std::vector<int>& bandlanczos<MATRIX, VS>::multiplicities() {
    return multiplicity;
  }
  
  template <class MATRIX, class VS>
    const std::vector<typename bandlanczos<MATRIX, VS>::magnitude_type>& bandlanczos<MATRIX, VS>::eigenvalues() {
    return eigenvals;
  }
   
  template <class MATRIX, class VS>
    const std::vector<typename bandlanczos<MATRIX, VS>::vector_type>& bandlanczos<MATRIX, VS>::eigenvectors() {
    return reigenvecs;
  }
  
  template <class MATRIX, class VS> template <class GEN, class ITER> 
    void bandlanczos<MATRIX, VS>::calculate_eigenvectors(ITER& iter, const GEN& gen) {
    magnitude_type rnorm;
    scalar_type cnorm;
    int j,k,k0,i;
    reigenvecs = std::vector<vector_type>(eigenvecs.size());
    for (i=0;i<reigenvecs.size();i++)
      reigenvecs[i] = new_vector(vecspace_);
    
    // set v_k starting vectors for k=1,...,p
    indexer index(p_);
    int pc=p_;
    v = std::vector<vector_type>(2*p_+1);
    for (i=0;i<(2*p_+1);i++)
      v[i] = new_vector(vecspace_);
    for (int x=0;x<p_;x++)
      ietl::generate(v[x],gen);
    
    // set I={}
    std::vector<int> I;
    for (j=0;(j<=iter.iterations() && pc>0);j++) {
      
      // compute |v_j|_2
      rnorm=ietl::two_norm(v[index.cnv(j)]);
      
      //  decide if v_j should be deflated
      if (rnorm < iter.def_tol()) {        
        // if (j-p_c) > 0 set I = I and {j-p_c}
        if ( (j-pc)>=0 ) {
          I.push_back(j-pc);
          index.deflate(j-pc);
        }
        
        // set p_c = p_c -1
         pc--;
         
         // if p_c = 0 set j=j-1 and stop
         if (pc==0)
           continue;
         else {             
           // for k=j,...,j+p_c-1  set  v_k = v_{k+1}
           for (k=j;k<j+pc;k++)
             std::swap(v[index.cnv(k)],v[index.cnv(k+1)]);
           
           // return to compute |v_j|_2
           j--;
           continue;
         }  // END OF (if (pc==0)
      } // END OF if (norm<dtol)
      
      // normalize v_j = v_j / t_{j,j-p_c}
      v[index.cnv(j)]/=rnorm;
      
      // for k=j+1,...,j+p_c-1  set v_k = v_k - v_j t_{j,k-p_c}
      for (k=j+1;k<j+pc;k++) {
        cnorm = ietl::dot(v[index.cnv(j)],v[index.cnv(k)]);
        v[index.cnv(k)] -=cnorm*v[index.cnv(j)];
      } // END OF for (k=j+1;k<j+pc;k++)
      
      // v_j now won't be modified anymore. It's a Krylov basis vector, thus can be used
      // to transform the eigenvectors in the Krylov basis back to the original basis
      for (i=0;i<eigenvecs.size();i++)
        reigenvecs[i] += eigenvecs[i][j]*v[index.cnv(j)];
      
      // compute v_{j+p_c} = A v_j
      ietl::mult(matrix_, v[index.cnv(j)], v[index.cnv(j+pc)]);
       
      // set k0 = max{1,j-p_c}. For k=k0,...,j-1 
      k0 = ( (j-pc)>0 ? j-pc : 0 );
      for (k=k0;k<j;k++) // set v_{j+p_c} = v_{j+pc} - v_k t_{k,j}
        v[index.cnv(j+pc)] -=tmat(k,j)*v[index.cnv(k)];
      
      // for k in I and {j} (in ascending order)
      std::sort(I.begin(),I.end());
      for (i=0;i<I.size()+1;i++) {           
        // set v_{j+p_c} = v_{j+p_c} - v_k t_{k,j}
        if (i==I.size()) 
          k=j; 
        else 
          k=I[i];
        cnorm=ietl::dot(v[index.cnv(k)],v[index.cnv(j+pc)]);
        v[index.cnv(j+pc)]-=cnorm*v[index.cnv(k)];
      } 
      
      index.next();
    } // END OF for (j=0;....)
    
    // normalize the eigenvectors which are now given in the original basis
    for (i=0;i<reigenvecs.size();i++)
      reigenvecs[i]/=ietl::two_norm(reigenvecs[i]);
  }
  
  template <class MATRIX, class VS> template <class GEN, class ITER>
    void bandlanczos<MATRIX, VS>::calculate_eigenvalues(ITER& iter, const GEN& gen) {
    magnitude_type rnorm;
    scalar_type cnorm;
    int j,k,k0,i;
    tmat.resize(iter.max_iter()+1, iter.max_iter()+1);
    
    // set v_k starting vectors for k=1,...,p
    indexer index(p_);
    for (i=0;i<(2*p_+1);i++)
      v[i] = new_vector(vecspace_);
    for (int x=0;x<p_;x++)
      ietl::generate(v[x],gen);
         
    // set p_c = p  and I = {}
    int pc=p_;
    std::vector<int> I;
    
    // start iteration until max_iterations reached or pc=0
    for (j=0;( (j<=iter.max_iter()) && (pc>0) );j++) {
      ++iter;
      
      // compute |v_j|_2
      rnorm=ietl::two_norm(v[index.cnv(j)]);
      
      // decide if v_j should be deflated
      if (rnorm < iter.def_tol()) {
        // if j-p_c > 0 set I=I and {j-p_c}
        if ((j-pc)>=0) {
          I.push_back(j-pc);
          index.deflate(j-pc);
        }
        
        // set p_c = p_c -1
        pc--;
        
        // if p_c = 0 set j=j-1 and stop
        if (pc==0)
          continue;
        else {                   
          // for k=j, j+1, ... , j+p_c-1  set  v_k = v_{k+1}
          for (k=j;k<j+pc;k++)
            std::swap(v[index.cnv(k)],v[index.cnv(k+1)]);
          
          // return to compute |v_j|_2
          j--;
          --iter;
          continue;
        }  // END OF (if (pc==0)
      } // END OF if (norm<dtol)
      
      // set t_{j, j-p_c} = |v_j|_2 and normalize v_j = v_j / t_{j, j-p_c}
      if (j-pc>=0)
        tmat(j,j-pc) = rnorm;
      v[index.cnv(j)] /= rnorm;
           
      // for k=j+1, ... , j+p_c -1
      for (k=j+1;k<j+pc;k++) {               
        // set t_{j, k-p_c} = v_j^\star v_k  and  v_k = v_k - v_j t_{j, k-p_c}
        cnorm = ietl::dot(v[index.cnv(j)],v[index.cnv(k)]);
        if (k-pc>=0)
          tmat(j,k-pc) = cnorm;
        v[index.cnv(k)] -= cnorm * v[index.cnv(j)];
      }
      
      // compute v_{j+p_c} = A v_j
      ietl::mult(matrix_, v[index.cnv(j)], v[index.cnv(j+pc)]);
      
      // set k_0 = max{1,j-p_c}. For k=k0,...,j-1
      k0 = ( (j-pc) > 0 ? j-pc : 0 );
      for (k=k0;k<j;k++) {               
              // set t_{k,j} = \bar{t_{j,k}} and v_{j+p_c} = v_{j+p_c} - v_k t_{k,j}
        tmat(k,j) = tmat(j,k);
        v[index.cnv(j+pc)] -= tmat(k,j) * v[index.cnv(k)];
      }
      
      // for k in I and {j} (in ascending order)
      std::sort(I.begin(),I.end());
      for (i=0;i<I.size()+1;i++) {             
        // set t_{k,j} = v_k^\star v_{j+p_c}  and  v_{j+p_c} = v_{j+p_c} - v_k t_{k,j}
        if (i==I.size()) k=j; else k=I[i];
        cnorm=ietl::dot(v[index.cnv(k)],v[index.cnv(j+pc)]);
        tmat(k,j) = cnorm;
        v[index.cnv(j+pc)] -= cnorm * v[index.cnv(k)];
      }
           
      // for k in I
      for (i=0;i<I.size();i++) {               
        // set t_{j,k} = \bar{t_{k,j}}
        k=I[i];
        tmat(j,k) = tmat(k,j);
      }
      
      // test for convergence
      if ( convergence_test(j+1, pc, iter.evs(), iter.dep_tol(), iter.ghost_tol(), iter.ghost_discarding(), iter.low()) == true )
        break;
      index.next();
    } // END OF for (j=0;....)
    calc_approx_eigenvalues(iter.iterations()+1, iter.evs(), iter.ghost_tol(), iter.ghost_discarding(), iter.low());
    // can be used for debugging, this returns ALL computed eigenvalues, even when not converged...
    // calc_approx_eigenvalues(iter.iterations()+1, -1, iter.ghost_tol(), iter.ghost_discarding(), iter.low());
  }
  
  template <class MATRIX, class VS>
    bool bandlanczos<MATRIX, VS>::convergence_test(const int dim, int pc, int evs, magnitude_type dep_tol, magnitude_type ghost_tol, bool ghost_discarding, bool low) {
    if (dim<evs)  // cannot have converged, because dim of tridiagonal matrix is smaller than #wanted ev 
      return false;
    int num_evs = calc_approx_eigenvalues(dim,-1, ghost_tol,ghost_discarding,low);
    if (num_evs < evs)        // cannot have converged, we have less eigenvalues than wanted
      return false;
    if ( multiplicity[evs-1]>1 )        // converged, since eigenvectors are multiple
      return true;
    magnitude_type norm=0;
    magnitude_type vnorm=0;
    for (int i=1;i<=evs;i++) {
      vnorm = 0;
      for (int j=1;j<=pc;j++)
        vnorm += std::abs(eigenvecs[evs-i][dim-j]);// * std::abs(eigenvecs[evs-i][dim-j]);
      norm = ( vnorm > norm ? vnorm : norm );
    }
    if ( (norm <= dep_tol) ) 
      return true;
    else
      return false;
  }
   
  template<class MATRIX, class VS>
    int bandlanczos<MATRIX, VS>::calc_approx_eigenvalues(const int& dim, int evs, magnitude_type ghost_tol, bool ghost_discarding, bool low) {
    eigenvals = std::vector< magnitude_type >();
    eigenvecs = std::vector< eigenvec_type >();
    multiplicity = std::vector< int >();
    eigenvec_type vec(dim);
    Tjpr.resize(dim,dim);
    for (int i=0;i<dim;i++) for (int j=0;j<=i;j++)
      Tjpr(j,i) = tmat(j,i);
    magnitude_type* w = new magnitude_type[dim];
    int info=0;
    scalar_type type=0;
    getev_(dim, w, info, type);
    int* mult = new int[dim];
    int* evpl = new int[dim];// Stores the index of w, where the eigenvals are (second ev is at place 7 -> evpl[7]=2)
    // bool ghost[dim];
    bool* ghost = new bool[dim];
    int num_evs;
    int discarded=0; // number of ghosts, is substracted from #eigenvalues of T to get #ev of the original matrix
    evpl[0] = 0;
    mult[0] = 1;
    ghost[0]= false;
    num_evs=1;
    for (int i=0;i<dim-1;i++) {
      if ( std::abs((w[i] - w[i+1]) ) > ghost_tol) {
        evpl[num_evs] = i+1;
        mult[num_evs] = 1;
        ghost[num_evs]= false;
        num_evs++;
      } 
      else
        mult[num_evs-1]++;
    }
    for (int i=0;i<num_evs;i++) {
      for (int j=0;j<dim;j++)
        if (low)
          vec[j] = Tjpr(j,evpl[i]);
        else
          vec[j] = Tjpr(j,evpl[num_evs-i-1]);  // CHANGE!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      eigenvecs.push_back(vec);
    }
    if (((num_evs>evs) && (ghost_discarding)) && (dim>1)) {
      Tjpr.resize(dim-1,dim-1);
      for (int i=1;i<dim;i++)
        for (int j=1;j<=i;j++)
          Tjpr(j-1,i-1) = tmat(j,i);
      magnitude_type* w2 = new magnitude_type[dim-1];
      getev_(dim-1,w2,info,type);
      for (int i=0;i<num_evs;i++) {
        if (mult[i]==1) {
          for (int j=0;j<dim-1;j++) {
            if ( std::abs( (w2[j] - w[evpl[i]]) ) < ghost_tol) {
              ghost[i]=true;
              if (low)
                eigenvecs.erase(eigenvecs.begin()+i-discarded,eigenvecs.begin()+i+1-discarded);
              else
                eigenvecs.erase(eigenvecs.begin()+num_evs-i-1,eigenvecs.begin()+num_evs-i); // CHANGE !!!
              discarded++;
            }
          }
        }
      }
      delete[] w2;
    }
    if (evs == -1)
      for (int i=0;i<num_evs;i++) {
        if (low){
          if (!ghost[i]) {
            eigenvals.push_back(w[evpl[i]]);
            multiplicity.push_back(mult[i]);
          }
        }
        else {
          if (!ghost[num_evs-i-1]) {                  
            eigenvals.push_back(w[evpl[num_evs-i-1]]); // CHANGE
            multiplicity.push_back(mult[num_evs-i-1]); // CHANGE
          }
        }
      }
    else {
      discarded = 0; // now used to check if we already stored the wanted eigenvalues
      for (int i=0;( (i<num_evs) && (discarded<evs) ); i++) {
        if ( low ) {
          if (!ghost[i]) {
            eigenvals.push_back(w[evpl[i]]);
            multiplicity.push_back(mult[i]);
            discarded++;
          }
           }
        else {
          if (!ghost[num_evs-i-1]) {
            eigenvals.push_back(w[evpl[num_evs-i-1]]);  // CHANGE
            multiplicity.push_back(mult[num_evs-i-1]);  // CHANGE
            discarded++;
          }
        }
      }
      eigenvecs.erase(eigenvecs.begin()+discarded,eigenvecs.end());
    } 
    delete[] w;
    delete[] mult;
    delete[] evpl;
    delete[] ghost;
    if (evs == -1)
      return num_evs-discarded;
    else
      return discarded;       
  }
      
  template<class MATRIX, class VS>
    void bandlanczos<MATRIX, VS>::getev_(const fortran_int_t dim, double w[], fortran_int_t info, double type) {
    char jobz='V';
    char uplo='U';
    const fortran_int_t in=dim;
    fortran_int_t lda=dim;
    const fortran_int_t lwork=4*in; 
    double* work = new double[lwork];
    LAPACK_DSYEV(jobz, uplo, in, Tjpr.data(), lda, w, work, lwork, info);
    delete[] work;
  }
  
  template<class MATRIX, class VS>
    void bandlanczos<MATRIX, VS>::getev_(const fortran_int_t dim, float w[], fortran_int_t info, float type){
    char jobz='V';
    char uplo='U';
    int in=dim;
    int lda=dim;
    int lwork=4*in;
    float* work = new float[lwork];
    LAPACK_SSYEV(jobz, uplo, in, Tjpr.data(), lda, w, work, lwork, info);
    delete[] work;
  }
  
  template<class MATRIX, class VS>
    void bandlanczos<MATRIX, VS>::getev_(const fortran_int_t dim, double w[], fortran_int_t info, std::complex<double> type) {
    char jobz='V';
    char uplo='U';
    fortran_int_t in=dim;
    fortran_int_t lda=dim;
    fortran_int_t lwork=4*in;
    std::complex<double> work[lwork];
    double rwork[lwork];
    LAPACK_ZHEEV(jobz, uplo, in, Tjpr.data(), lda, w, work, lwork, rwork, info);
  }
      
  template<class MATRIX, class VS>
    void bandlanczos<MATRIX, VS>::getev_(const fortran_int_t dim, float w[], fortran_int_t info, std::complex<float> type) {
    char jobz='V';
    char uplo='U';
    fortran_int_t in=dim;
    fortran_int_t lda=dim;
    fortran_int_t lwork=4*in;
    std::complex<float> work[lwork];
    float rwork[lwork];
    LAPACK_CHEEV(jobz, uplo, in, Tjpr.data(), lda, w, work, lwork, rwork, info);
  }  
  // E N D   O F   C L A S S   B A N D L A N C Z O S /////////////////////////
}
#endif
