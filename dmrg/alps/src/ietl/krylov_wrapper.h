/*****************************************************************************
*
* ALPS Project: Algorithms and Libraries for Physics Simulations
*
* ALPS Libraries
*
* Copyright (C) 2002 by Matthias Troyer <troyer@comp-phys.org>
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

#ifndef KRYLOV_WRAPPER_H
#define KRYLOV_WRAPPER_H

#include <itl/krylov/cg.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >  

class cg_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;
  
 public:
  cg_wrapper(const Preconditioner M, Iteration iter)
    : M_(M), iter_(iter), N_(0)
  {
  }

template < class scalar_type >  
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    { 
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b 
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration cg_iter = iter_;
    int cg_error_code = cg(m2, x, b, M_, cg_iter);
    if (cg_error_code != 0)
    {
      throw std::runtime_error("ERROR in CG_WRAPPER:    cg can't solve system");
    }
  }
};
}


#include <itl/krylov/cgs.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >

class cgs_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;

 public:
  cgs_wrapper(const Preconditioner M, Iteration iter)
    : M_(M), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int cg_error_code = cgs(A, x, b, M_, iter);
    if (cg_error_code != 0)
    {
      throw std::runtime_error("ERROR in CGS_WRAPPER:    cgs can't solve system");
    }
  }
};
}


#include <itl/krylov/bicg.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >  

class bicg_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;
 
 public:
  bicg_wrapper(const Preconditioner M, Iteration iter)
    : M_(M), iter_(iter), N_(0) 
  {
  }

  template < class scalar_type >  
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    { 
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b 
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int cg_error_code = bicg(A, x, b, M_, iter);
    if (cg_error_code != 0)
    {
      throw std::runtime_error("ERROR in BICG_WRAPPER:    bicg can't solve system");
    }
  }
};
}

/*
#include <itl/krylov/gmres.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration, class Basis >
class gmres_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Basis KS_;
  int restart_;
  Matrix unit_;
  int N_;

 public:
  gmres_wrapper(const Preconditioner M, Iteration iter, int restart, Basis KS)
    : M_(M), iter_(iter), KS_(KS), restart_(restart), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code =gmres(A, x, b, p, restart, iter, KS);
    if (error_code != 0)
    {
      throw std::runtime_error("ERROR in GMRES_WRAPPER:    bicg can't solve system");
    }
  }
};
}
*/

#include <itl/krylov/bicgstab.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >

class bicgstab_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;

 public:
  bicgstab_wrapper(const Preconditioner M, Iteration iter)
    : M_(M), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code = bicgstab(A, x, b, M_, iter);
    if (error_code != 0)
    {
      throw std::runtime_error("ERROR in BICGSTAB_WRAPPER:    bicgstab can't solve system");
    }
  }
};
}


#include <itl/krylov/qmr.h>

namespace ietl {
template <class Matrix, class VectorX, class VectorB, class Precond1, class Precond2, class Iteration>

class qmr_wrapper
{
 private:
  Precond1 M1_;
  Precond2 M2_;
  Iteration iter_;
  Matrix unit_;
  int N_;

 public:
  qmr_wrapper(const Precond1 M1, const Precond2 M2, Iteration iter)
    : M1_(M1), M2_(M2), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int qmr_error_code = qmr(A, x, b, M1_, M2_, iter);
    if (qmr_error_code != 0)
    {
      switch (qmr_error_code)
      {
      case 1:
          throw std::runtime_error("ERROR in QMR_WRAPPER: no convergence after maximum iterations");
      case 2:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in rho");
      case 3:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in beta");
      case 4:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in gamma");
      case 5:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in delta");
      case 6:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in ep");
      case 7:
          throw std::runtime_error("ERROR in QMR_WRAPPER: breakdown in x");
      default:
        throw std::runtime_error("ERROR in QMR_WRAPPER");
      }
    }
  }
};
}


#include <itl/krylov/tfqmr.h>

namespace ietl {
template <class Matrix, class VectorX, class VectorB, class Precond1, class Precond2, class Iteration>

class tfqmr_wrapper
{
 private:
  Precond1 M1_;
  Precond2 M2_;
  Iteration iter_;
  Matrix unit_;
  int N_;

 public:
  tfqmr_wrapper(const Precond1 M1, const Precond2 M2, Iteration iter)
    : M1_(M1), M2_(M2), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code = tfqmr(A, x, b, M1_, M2_, iter);
    if (error_code != 0)
    {
      switch (error_code)
      {
      case 1:
          throw std::runtime_error("ERROR in TFQMR_WRAPPER: no convergence after maximum iterations");
      case 2:
          throw std::runtime_error("ERROR in TFQMR_WRAPPER: breakdown in tau");
      case 3:
          throw std::runtime_error("ERROR in TFQMR_WRAPPER: breakdown in alpha");
      case 4:
          throw std::runtime_error("ERROR in TFQMR_WRAPPER: breakdown in gamma");
      case 5:
          throw std::runtime_error("ERROR in TFQMR_WRAPPER: breakdown in rho");
      default:
        throw std::runtime_error("ERROR in TFQMR_WRAPPER");
      }
    }
  }
};
}


#include <itl/krylov/gcr.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >

class gcr_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  int restart_;
  Matrix unit_;
  int N_;

 public:
  gcr_wrapper(const Preconditioner M, Iteration iter, int restart)
    : M_(M), iter_(iter), restart_(restart), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code =gcr(A, x,  b, M_, restart_, iter);
    if (error_code != 0)
    {
      throw std::runtime_error("ERROR in GCR_WRAPPER:    gcr can't solve system");
    }
  }
};
}


#include <itl/krylov/cheby.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >

class cheby_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;

  typename VectorB::value_type eigmax_;
  typename VectorB::value_type eigmin_;

 public:
  cheby_wrapper(const Preconditioner M, Iteration iter, typename VectorB::value_type eigmin, typename VectorB::value_type eigmax)
    : M_(M), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code = cheby(A, x, b, M_, iter);
    if (error_code != 0)
    {
      throw std::runtime_error("ERROR in CHEBY_WRAPPER:    cheby can't solve system");
    }
  }
};
}


#include <itl/krylov/richardson.h>

namespace ietl {
template < class Matrix, class VectorX, class VectorB, class Preconditioner, class Iteration >

class richardson_wrapper
{
 private:
  Preconditioner M_;
  Iteration iter_;
  Matrix unit_;
  int N_;

 public:
  richardson_wrapper(const Preconditioner M, Iteration iter)
    : M_(M), iter_(iter), N_(0)
  {
  }

  template < class scalar_type >
  void operator()(const Matrix& A, scalar_type s, VectorX& x, const VectorB& b) throw (std::runtime_error)
  {
    if (N_ == 0)
    {
      N_ = A.ncols();
      // generate identity-Matrix
      Matrix unit(N_, N_);
      for(int i = 0; i < N_; i++)
      {
        unit(i,i) = 1;
      }
      unit_ = unit;
    }

    // working variables
    Matrix m1(N_, N_);
    Matrix m2(N_, N_);

    // calc (A-s*eye(N))/x = b
    itl::copy(unit_, m1);
    itl::scale(m1, -s);
    itl::copy(A, m2);
    itl::add(m1, m2);

    Iteration iter = iter_;
    int error_code = richardson(A, x, b, M_, iter);
    if (error_code != 0)
    {
      throw std::runtime_error("ERROR in RICHARDSON_WRAPPER:    richardson can't solve system");
    }
  }
};
}


#endif

