/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_LANCZOSEVOLVER_H
#define MAQUIS_DMRG_LANCZOSEVOLVER_H

#ifdef DMRG_TD

#include <vector>
//#include <Eigen/Core>
//#include <unsupported/Eigen/MatrixFunctions>
#include "TimeEvolutionAlgorithm.h"

enum class TimeStepDistributor { ExponentialMidpoint, FourthOrderMagnus, Standard };

/**
 * @brief Trait class that defines the intermediate time-steps and coefficients 
 * for all available Lanczos-type propagators.
 */

template<TimeStepDistributor TimeStepDistributorType>
struct TimeStepTraits {};

template<>
struct TimeStepTraits<TimeStepDistributor::ExponentialMidpoint> {
  static constexpr int numberOfExponentials = 1;
  static constexpr int numberOfFactorsPerExponential = 1;
  using InternalValueType = std::array<std::pair<double, double>, 1>;
  using FactorsType = std::array<InternalValueType, numberOfExponentials>;
  static constexpr FactorsType factorsAndSteps = {{std::make_pair(1., 0.5)}};
};

template<>
struct TimeStepTraits<TimeStepDistributor::FourthOrderMagnus> {
  static constexpr int numberOfExponentials = 2;
  static constexpr int numberOfFactorsPerExponential = 2;
  using InternalValueType = std::array<std::pair<double, double>, 2>;
  using FactorsType = std::array<InternalValueType, numberOfExponentials>;
  static constexpr FactorsType factorsAndSteps 
    = { InternalValueType({std::make_pair((3.-2.*sqrt(3.))/12., 0.5+sqrt(3.)/6.), std::make_pair((3.+2.*sqrt(3.))/12., 0.5-sqrt(3.)/6.)}),
        InternalValueType({std::make_pair((3.+2.*sqrt(3.))/12., 0.5+sqrt(3.)/6.), std::make_pair((3.-2.*sqrt(3.))/12., 0.5-sqrt(3.)/6.)}) };
};

template<>
struct TimeStepTraits<TimeStepDistributor::Standard> {
  static constexpr int numberOfExponentials = 1;
  static constexpr int numberOfFactorsPerExponential = 1;
  using InternalValueType = std::array<std::pair<double, double>, 1>;
  using FactorsType = std::array<InternalValueType, numberOfExponentials>;
  static constexpr FactorsType factorsAndSteps = {{std::make_pair(1., 1.)}};
};

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
class LanczosEvolver : public TimeEvolutionAlgorithm<Matrix, SymmGroup> {

  /** Types definition */
  using complex_type = std::complex<double>;
  using base = typename TimeEvolutionAlgorithm<Matrix, SymmGroup>::TimeEvolutionAlgorithm;
  using scalar_type = typename MPSTensor<Matrix, SymmGroup>::scalar_type;
  using time_type = typename base::time_type;
  //using matrix_complex = Eigen::Matrix< complex_type, Eigen::Dynamic, Eigen::Dynamic >;
  //using vector_complex = Eigen::Matrix< complex_type, Eigen::Dynamic, 1 >;
  using matrix_complex = alps::numeric::matrix< complex_type >;
  using vector_complex = std::vector< complex_type >;
  using FactorsType = typename TimeStepTraits<TimeStepDistributorClass>::FactorsType;

  /** Types inheritance */
  using base::has_td_part_;
  using base::is_imag_;
  using base::time_step_;
  using base::apply_hamiltonian;

 public:

  /* Class constructor */
  LanczosEvolver(time_type time_step, bool has_td, bool is_imag, double threshold, std::size_t max_iter) 
    : base(time_step, has_td, is_imag), threshold_(threshold), max_iter_(max_iter) {}

  /* Time evolution method */
  void evolve(SiteProblem<Matrix, SymmGroup> const& site_problem, MPSTensor<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const override final {
    evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

  void evolve(ZeroSiteProblem<Matrix, SymmGroup> const& site_problem, block_matrix<Matrix, SymmGroup>& matrix,
              bool is_forward, time_type time_current, time_type time_step) const override final {
    evolve_kernel(site_problem, matrix, is_forward, time_current, time_step);
  }

 private:

  /* Kernel for the time evolution part */
  template<class SiteProblem, class MatrixType>
  void evolve_kernel(SiteProblem const& site_problem, MatrixType& matrix, bool is_forward, time_type time_current,
                     time_type time_step) const;

  /* Private method interfacing to Eigen matrix exponential calculator */
  template<class MatrixType, class VectorType>
  void apply_exponential(MatrixType& hamiltonian_matrix, VectorType& ret, size_t local_dim_, bool isForward) const;
  template<class SiteProblem, class MatrixType>
  MatrixType applyOperator(const MatrixType& inputVec, const SiteProblem& site_problem, int idExp, time_type time_current) const;

  /* Real --> Complex conversion routines */
  template< class ArgType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  complex_type initial_convert(const ArgType& input) const { return std::complex<double>(input, 0.) ; };
  template< class ArgType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  complex_type initial_convert(const ArgType& input) const { return input ; };
  template< class ArgType, typename std::enable_if< std::is_same<double, ArgType >::value>::type * = nullptr >
  ArgType final_convert(const complex_type& input) const { return std::real(input) ; };
  template< class ArgType, typename std::enable_if< std::is_same< typename std::complex<double>, ArgType >::value>::type * = nullptr >
  ArgType final_convert(const complex_type& input) const { return input ; };

  /* Methods to print the results of the Lanczos algorithm. */
  void print_header() const {
    print_line() ;
    std::cout << "  ITERATION  |   ERROR " << std::endl ;
    print_line() ;
  }
  void print_line() const {
    std::cout << "+------------+------------+" << std::endl ;
  }
  void print_data(std::size_t n_iter, time_type error) const {
    char buf[100] ;
    int n = sprintf(buf, "   %2d        |  %1.4E ", static_cast<int>(n_iter), error) ;
    std::cout << buf << std::endl ;
  };

  /* Class members */
  std::size_t max_iter_;
  double threshold_;
  static const FactorsType factorsAndSteps;
  static const int numberOfExponentials;
  static const int numberOfFactorsPerExponential;
};

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
constexpr int LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::numberOfExponentials
  = TimeStepTraits<TimeStepDistributorClass>::numberOfExponentials;

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
constexpr int LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::numberOfFactorsPerExponential
  = TimeStepTraits<TimeStepDistributorClass>::numberOfFactorsPerExponential;

template<class Matrix, class SymmGroup, TimeStepDistributor TimeStepDistributorClass>
constexpr typename LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::FactorsType
  LanczosEvolver<Matrix, SymmGroup, TimeStepDistributorClass>::factorsAndSteps = TimeStepTraits<TimeStepDistributorClass>::factorsAndSteps;

#include "LanczosEvolver.cpp"

#endif // DMRG_TD

#endif // MAQUIS_DMRG_RKEVOLVER_H
