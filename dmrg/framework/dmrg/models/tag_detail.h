/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2013-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *****************************************************************************/

#ifndef MAQUIS_DMRG_MODELS_TAG_DETAIL_H
#define MAQUIS_DMRG_MODELS_TAG_DETAIL_H

#include <alps/numeric/isnan.hpp>
#include <alps/numeric/isinf.hpp>

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
    bool is_uniform(block_matrix<Matrix, SymmGroup> const& op)
    {
        typename Matrix::value_type invscale;
        #ifdef AMBIENT
        {
            locale_shared l;
            for(int i = 0; i < op.n_blocks(); ++i) ambient::migrate(op[i]);
            ambient::sync();
        }
        #endif
        
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
                typename Matrix::value_type normalized = std::abs(m(i,j) * invscale);
                // if not 1 and not 0
                if (std::abs(normalized-1.0) > 1e-15 && normalized > 1e-15)
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
        #ifdef AMBIENT
        {
            locale_shared l;
            for(int i = 0; i < reference.n_blocks(); ++i) ambient::migrate(reference[i]);
            for(int i = 0; i < sample.n_blocks(); ++i) ambient::migrate(sample[i]);
            ambient::sync();
        }
        #endif
        if (reference.left_basis() != sample.left_basis() || reference.right_basis() != sample.right_basis())
            return std::make_pair(false, 0.);

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
