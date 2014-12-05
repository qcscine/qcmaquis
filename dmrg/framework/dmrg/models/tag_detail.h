/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 * 
 * This software is part of the ALPS Applications, published under the ALPS
 * Application License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 * 
 * You should have received a copy of the ALPS Application License along with
 * the ALPS Applications; see the file LICENSE.txt. If not, the license is also
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

#ifndef MAQUIS_DMRG_MODELS_TAG_DETAIL_H
#define MAQUIS_DMRG_MODELS_TAG_DETAIL_H

#include <alps/numeric/isnan.hpp>
#include <alps/numeric/isinf.hpp>
#include <alps/numeric/is_nonzero.hpp>

namespace tag_detail {

    typedef unsigned tag_type;

    enum operator_kind { bosonic, fermionic };

    struct pair_cmp
    {
        bool operator()(std::pair<tag_type, tag_type> const & i,
                        std::pair<tag_type, tag_type> const & j) const
        {
            if (i.first < j.first)
                return true;
            else if (i.first > j.first)
                return false;
            else
                return i.second < j.second;
        }
    };

    template <class Matrix, class SymmGroup>
    void remove_empty_blocks(block_matrix<Matrix, SymmGroup> & op)
    {
        {
            parallel::guard::serial guard;
            storage::migrate(op);
        }

        for (typename Matrix::size_type b=0; b < op.n_blocks(); ++b)
        {
            bool only_zero = true;
            const Matrix& m = op[b];
            for (int i = 0; i < num_rows(m); i++)
               for(int j = 0; j < num_cols(m); j++)
            {
                if (alps::numeric::is_nonzero(m(i,j))) {
                    only_zero = false;
                    break;
                }
            }
            if (only_zero) {
                op.remove_block(b);
                --b;
            }
        }
    }

    template <class Matrix, class SymmGroup>
    bool is_uniform(block_matrix<Matrix, SymmGroup> const& op)
    {
        if (op.n_blocks() == 0)
            return true;

        typename Matrix::value_type invscale;
        {
            parallel::guard::serial guard;
            storage::migrate(op);
        }
        
        // determine scale of matrices
        const Matrix& m = op[0];
        for (int i = 0; i < num_rows(m); i++)
           for(int j = 0; j < num_cols(m); j++)
            if (std::abs(m(i,j)) > 1.e-20) {
                invscale = 1./m(i,j);
                break;
            }

        for (typename Matrix::size_type b=0; b < op.n_blocks(); ++b)
        {
            const Matrix& m = op[b];
            for (int i = 0; i < num_rows(m); i++)
               for(int j = 0; j < num_cols(m); j++)
            {
                typename maquis::traits::real_type<typename Matrix::value_type>::type normalized = std::abs(m(i,j) * invscale);
                // if not 1 and not 0
                if (std::abs(normalized-1.0) > 1e-15 && std::abs(normalized) > 1e-15)
                    return false;
            }
        }

        return true;
    }

    template <class T>
    bool num_check(T x) {
        if (alps::numeric::isnan(x) || alps::numeric::isinf(x))
            throw std::runtime_error("NaN / INF numeric Error occured while comparing operator scales\n");
        return true;
    }

    inline bool num_check(std::complex<double> x) { return true; }

    template <class Matrix, class SymmGroup>
    std::pair<bool, typename Matrix::value_type>
    equal(block_matrix<Matrix, SymmGroup> const& reference,
          block_matrix<Matrix, SymmGroup> const& sample)
    {
        {
            parallel::guard::serial guard;
            storage::migrate(reference);
            storage::migrate(sample);
        }

        if (!shape_equal(reference, sample))
            return std::make_pair(false, 0.);

        if (sample.n_blocks() == 0)
            return std::make_pair(true, 1.0);

        typename Matrix::value_type invscale1, invscale2;
     
        // determine scale of matrices
        const Matrix& m1 = reference[0];
        for (int i = 0; i < num_rows(m1); i++)
           for(int j = 0; j < num_cols(m1); j++)
        {
            if (std::abs(m1(i,j)) > 1.e-50) {
                invscale1 = 1./m1(i,j);
                break;
            }
            if(i == (num_rows(m1)-1) && j == (num_cols(m1)-1)){ return std::make_pair(false, 0.); }
        }

        const Matrix& m2 = sample[0];
        for (int i = 0; i < num_rows(m2); i++)
           for(int j = 0; j < num_cols(m2); j++)
        {
            if (std::abs(m2(i,j)) > 1.e-50) {
                invscale2 = 1./m2(i,j);
                break;
            }
            if(i == (num_rows(m2)-1) && j == (num_cols(m2)-1)){ return std::make_pair(false, 0.); }
        }

        // Check all blocks for equality modulo scale factor
        for (typename Matrix::size_type b=0; b < reference.n_blocks(); ++b)
        {
            const Matrix& mb1 = reference[b];
            const Matrix& mb2 = sample[b];
            for (int i = 0; i < num_rows(mb1); i++)
               for(int j = 0; j < num_cols(mb1); j++)
            {
                typename Matrix::value_type t1 = mb1(i,j) * invscale1, t2 = mb2(i,j) * invscale2;
                if (std::abs(t1 - t2) > 1e-12)
                    return std::make_pair(false, 0.);
            }
        }

        typename Matrix::value_type scale = invscale1 / invscale2;

#ifndef NDEBUG
        try { num_check(invscale1); }
        catch (std::exception e) { maquis::cout << "invscale1 numcheck failed\n"; exit(1);}
        try { num_check(invscale2); }
        catch (std::exception e) { maquis::cout << "invscale2 numcheck failed\n"; exit(1);}
        try { num_check(scale); }
        catch (std::exception e) { maquis::cout << "scale numcheck failed\n"; exit(1);}
#endif

        return std::make_pair(true, scale);
    }
}

#endif
