
#ifndef BLOCK_MATRIX_ITERABLE
#define BLOCK_MATRIX_ITERABLE

#include <alps/numeric/real.hpp>
#include "utils/types.h"

namespace detail {
    
    template <class T>
    typename utils::real_type<T>::type gather_real_pred(T const & val)
    {
        assert( check_real(val) );
        assert( alps::numeric::real(val) > -1e-10 );
        return alps::numeric::real(val);
    }
    
    template <class DiagMatrix, class SymmGroup>
    struct iteretable_diag_impl {
        static void solver_truncate_impl_zero(block_matrix<DiagMatrix, SymmGroup> const & evals, size_t Mmax, double cutoff, size_t* keeps, double & truncated_weight, double & smallest_ev)
        {
            size_t length = 0;
            for(std::size_t k = 0; k < evals.n_blocks(); ++k){ 
                length += num_rows(evals[k]);
            }

            std::vector<typename utils::real_type<typename DiagMatrix::value_type>::type > allevals(length);

            std::size_t position = 0;
            for(std::size_t k = 0; k < evals.n_blocks(); ++k){
                std::transform(evals[k].elements().first, evals[k].elements().second, allevals.begin()+position, gather_real_pred<typename DiagMatrix::value_type>);
                position += num_rows(evals[k]);
            }

            assert( allevals.size() > 0 );
            std::sort(allevals.begin(), allevals.end());
            std::reverse(allevals.begin(), allevals.end());

            double evalscut = cutoff * allevals[0];

            if (allevals.size() > Mmax)
                evalscut = std::max(evalscut, allevals[Mmax]);
            smallest_ev = evalscut / allevals[0];
           
            truncated_weight = std::accumulate(std::find_if(allevals.begin(), allevals.end(), boost::lambda::_1 < evalscut), allevals.end(), 0.0);
            truncated_weight /= std::accumulate(allevals.begin(), allevals.end(), 0.0);
           
            for(std::size_t k = 0; k < evals.n_blocks(); ++k){
                std::vector<typename utils::real_type<typename DiagMatrix::value_type>::type> evals_k;
                for (typename DiagMatrix::const_element_iterator it = evals[k].elements().first; it != evals[k].elements().second; ++it)
                    evals_k.push_back(alps::numeric::real(*it));
                keeps[k] = std::find_if(evals_k.begin(), evals_k.end(), boost::lambda::_1 < evalscut)-evals_k.begin();
            }
        }
    };
}
#endif
