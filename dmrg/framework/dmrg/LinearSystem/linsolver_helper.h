/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef LINSOLVER_HELPER_H
#define LINSOLVER_HELPER_H

#include <tuple>

namespace LinSolverHelper {

/** @brief Small helper class representing a givens rotation */
template<class ScalarType>
class Givens {
public:
    
    /** @brief Default constructor */        
    explicit Givens() : activated_(false), c_(0.) {};

    /** @brief Constructor from a pair */
    Givens(ScalarType x0, ScalarType x1) : activated_(true) {
        std::tie(c_, s_) = calculateRotation(x0, x1);
        r_ = c_*x0 + s_*x1;
    }

    std::pair<ScalarType, ScalarType> apply(ScalarType x0, ScalarType x1) {
        return std::make_pair(c_*x0+s_*x1, -localConj(s_)*x0+c_*x1);
    }
    
    std::pair<double, std::complex<double>> calculateRotation(std::complex<double> ca, std::complex<double> cb) {
       if (std::abs(ca) < 1.0E-16) {
            c_ = 0.0;
            s_ = std::complex<double>(1., 0.); 
       }
       else {
          auto scale = std::abs(ca) + std::abs(cb);
          auto norm = scale*std::sqrt(std::pow(std::abs(ca/std::complex<double>(scale, 0.0)), 2)
                                    + std::pow(std::abs(cb/std::complex<double>(scale, 0.0)), 2));
          auto alpha = ca/std::abs(ca);
          c_ = std::abs(ca)/norm;
          s_ = alpha*localConj(cb)/norm;
       }
       return std::make_pair(c_, s_);
    }
    
    std::pair<double, double> calculateRotation(double da, double db) {
        auto roe = db;
        if (std::abs(da) > std::abs(db))
            roe = da;
        auto scale = std::abs(da) + std::abs(db);
        if (scale < 1.0E-16) {
            c_ = 1.0;
            s_ = 0.0;
            r_ = 0.0;
        }
        else { 
          r_ = scale*std::sqrt(std::pow(da/scale, 2) + std::pow(db/scale, 2));
          if (roe < 0.)
            r_ *= -1.;
          c_ = da/r_;
          s_ = db/r_;
        }
        return std::make_pair(c_, s_);
    }

    /** @brief Overloading of user-defined conjugate function */
    static inline double localConj(double a) {
        return a;
    }
    
    /** @brief Conjugate function for complex numbers */
    static inline std::complex<double> localConj(std::complex<double> a) {
        return std::conj(a);
    }

    bool isActivated() const {
        return activated_;
    }

    auto getR() const {
        return r_;
    }


private:
    double c_;
    ScalarType s_, r_;
    bool activated_;
};

} // LinsolverHelper        

#endif
