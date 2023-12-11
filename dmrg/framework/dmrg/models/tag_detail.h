/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Department of Chemistry and Applied Biosciences, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef MAQUIS_DMRG_MODELS_TAG_DETAIL_H
#define MAQUIS_DMRG_MODELS_TAG_DETAIL_H

#include <alps/numeric/isnan.hpp>
#include <alps/numeric/isinf.hpp>
#include <alps/numeric/is_nonzero.hpp>

#include "dmrg/utils/storage.h"

namespace tag_detail {

    using tag_type = unsigned int;

    enum operator_kind { bosonic, fermionic };

    template <class BlockMatrix>
    void remove_empty_blocks(BlockMatrix & op)
    {
        {
            parallel::guard::serial guard;
            storage::migrate(op);
        }

        for (typename BlockMatrix::size_type b=0; b < op.n_blocks(); ++b) {
          bool only_zero = true;
          for (int i = 0; i < num_rows(op[b]); i++) {
            for(int j = 0; j < num_cols(op[b]); j++) {
              if (alps::numeric::is_nonzero(op[b](i,j))) {
                only_zero = false;
                break;
              }
            }
          }
          if (only_zero) {
            op.remove_block(b);
            --b;
          }
        }
    }

    template <class T>
    bool num_check(T x) {
        if (alps::numeric::isnan(x) || alps::numeric::isinf(x)) {
            throw std::runtime_error("NaN / INF numeric Error occured while comparing operator scales\n");
        }
        return true;
    }

    inline bool num_check(std::complex<double> x) { return true; }

    
    /**
     * @brief Checks if two BlockMatrices are equal up to a scale factor
     *
     * @return Returns a pair containing either false and 0.
     *         or true and the scale factor
     */
    template <class BlockMatrix>
    std::pair<bool, typename BlockMatrix::matrix_type::value_type>
    is_scaled(BlockMatrix const& reference, BlockMatrix const& sample)
    {
        using Matrix = typename BlockMatrix::matrix_type;
        using value_type = typename Matrix::value_type;
 
        {
            parallel::guard::serial guard;
            storage::migrate(reference);
            storage::migrate(sample);
        }
        // if shapes are different, so are the matrices
        if (!shape_equal(reference, sample)) {
            return std::make_pair(false, 0.);
        }

        // if both have 0 blocks, they are empty and therefore equal
        if (sample.n_blocks() == 0) {
            return std::make_pair(true, 1.0);
        }

        // determine scale factor of the reference matrix
        const Matrix& m1 = reference[0];
        bool foundNonZero = false;
        // loop over all matrix elements until nonzero entry is found
        int i = 0, j = 0;
        for (; i < num_rows(m1); i++) {
            for(; j < num_cols(m1); j++) {
                if (std::abs(m1(i,j)) > 1.e-50) {
                    foundNonZero = true;
                    break;
                }
            }
            if (foundNonZero) break;
        }

        // Check whether the sample matrix also has a nonzero entry there
        const Matrix& m2 = sample[0];
        if (foundNonZero && std::abs(m2(i,j)) < 1.e-50) {
            return std::make_pair(false, 0.);
        }

        value_type invscale1 = value_type(1.);
        value_type invscale2 = value_type(1.);

        // If there is a nonzero entry, calculate the inverse scaling factors
        if (foundNonZero) {
            invscale1 = value_type(1.)/m1(i,j);
            invscale2 = value_type(1.)/m2(i,j);
        }


        // Check all blocks for equality modulo scale factor
        for (typename Matrix::size_type b=0; b < reference.n_blocks(); ++b) {
            const Matrix& mb1 = reference[b];
            const Matrix& mb2 = sample[b];
            for (int i = 0; i < num_rows(mb1); i++) {
              for(int j = 0; j < num_cols(mb1); j++) {
                typename Matrix::value_type t1 = mb1(i,j) * invscale1;
                typename Matrix::value_type t2 = mb2(i,j) * invscale2;
                if (std::abs(t1 - t2) > 1e-12) { 
                  return std::make_pair(false, 0.);
                }
              }
            }
        }

        typename Matrix::value_type scale = invscale1 / invscale2;

#ifndef NDEBUG
        try { num_check(invscale1); }
        catch (const std::exception& e) { maquis::cout << "invscale1 numcheck failed\n"; exit(1);}
        try { num_check(invscale2); }
        catch (const std::exception& e) { maquis::cout << "invscale2 numcheck failed\n"; exit(1);}
        try { num_check(scale); }
        catch (const std::exception& e) { maquis::cout << "scale numcheck failed\n"; exit(1);}
#endif

        return std::make_pair(true, scale);
    }
}

#endif
