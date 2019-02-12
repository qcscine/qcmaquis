/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2018         Leon Freitag <lefreita@ethz.ch>
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
#ifndef MULTI_CANONIZE_HPP
#define MULTI_CANONIZE_HPP

/*
    Here we implement a simultaneous canonicalisation of MPS that have been
    optimised simultaneously
*/


// multi_normalize_left()
// Split a vector of MPSTensors into a single left-normalised MPSTensor and a vector of MPSTensor
// which can be premultiplied by MPSTensors at the next site to obtain a mixed-canonical form of
// several MPS with the same basis
//
// inspired by the SVD of sa_alg = -3 by Alberto
//
// Input: std::vector<MPS<Matrix, SymmGroup> > mps_vec
// Output: RHS as the return value (to mimick MPSTensor<>::normalize_left()! )
// mps_vec will be normalised and all MPSTensors will be identical!
// TODO: This is implemented using std::reference_wrappers, but is not generic at all. Make it more generic!
template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > multi_normalize_left(std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > >& mps_vec, DecompMethod method)
{

    block_matrix<Matrix, SymmGroup> bm_overall;

    // make sure the normalisation of all MPSTensors is the same
    Indicator normalization = mps_vec[0].get().cur_normalization;
    for (auto&& mpsr : mps_vec)
        assert(mpsr.get().cur_normalization == normalization);

    // TODO: Add checks to make sure we've optimised MPS with the same boundaries.
    if (normalization == Unorm || normalization == Rnorm) {

        std::vector<block_matrix<Matrix, SymmGroup> > ret;
        for (auto&& mpsref : mps_vec)
        {
            MPSTensor<Matrix, SymmGroup> & mps = mpsref.get();
            // left-pair all MPSTensors
            mps.make_left_paired();
            // Adds the block to bm_overall
            for (std::size_t k = 0; k < mps.data().n_blocks(); ++k)
                bm_overall.add_block_to_column(mps.data(),
                                            mps.data().basis().left_charge(k),
                                            mps.data().basis().right_charge(k));
        }
        if (method == QR)
        {
            block_matrix<Matrix, SymmGroup> Q, R;
            // QR decompose bm_overall
            qr(bm_overall, Q, R);

            ret.reserve(mps_vec.size());

            for (auto&& mpsref : mps_vec)
            {
                block_matrix<Matrix, SymmGroup> tmp;
                MPSTensor<Matrix, SymmGroup> & mps = mpsref.get();
                // prepare the RHS
                // Premultiply Q^\dagger M on the std::vector<MPS<Matrix, SymmGroup> > elements
                // mps.make_right_paired(); //really?
                gemm(adjoint(Q), mps.data(), tmp);
                ret.push_back(tmp);

                // Now that we're ready with RHS, replace all MPSTensors with Q
                mps.data() = Q;
                mps.right_i = mps.data().right_basis();
                assert(mps.right_i == R.left_basis());
                mps.cur_normalization = Lnorm;
            }

            return ret;
        }
        else // if (method == SVD)
        {
            throw std::runtime_error("Multi_normalize_left with SVD NYI. Too lazy!!!");
        }

    }
    // assuming all elements in mps_vec have the same basis!!!
    // TODO: Add appropriate assertions
    return std::vector<block_matrix<Matrix, SymmGroup> >(mps_vec.size(), identity_matrix<block_matrix<Matrix, SymmGroup> >(mps_vec[0].get().data().right_basis()));

}

// multi_normalize_right()
// Same as multi_normalize_left(), but for right normalisation
template<class Matrix, class SymmGroup>
std::vector<block_matrix<Matrix, SymmGroup> > multi_normalize_right(std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > >& mps_vec, DecompMethod method)
{
    block_matrix<Matrix, SymmGroup> bm_overall;

    // make sure the normalisation of all MPSTensors is the same
    Indicator normalization = mps_vec[0].get().cur_normalization;
    for (auto&& mpsref : mps_vec)
        assert(mpsref.get().cur_normalization == normalization);

    if (normalization == Unorm || normalization == Lnorm) {

        std::vector<block_matrix<Matrix, SymmGroup> > ret;
        for (auto&& mpsref : mps_vec)
        {
            MPSTensor<Matrix, SymmGroup> & mps = mpsref.get();
            // right-pair all MPSTensors
            mps.make_right_paired();
            // Adds the block to bm_overall
            for (std::size_t k = 0; k < mps.data().n_blocks(); ++k)
                bm_overall.add_block_to_row(mps.data(),
                                            mps.data().basis().left_charge(k),
                                            mps.data().basis().right_charge(k));
        }
        if (method == QR) //enum QR but LQ decomposition
        {
            block_matrix<Matrix, SymmGroup> L, Q;
            // LQ decompose bm_overall
            lq(bm_overall, L, Q);

            ret.reserve(mps_vec.size());

            for (auto&& mpsref : mps_vec)
            {
                block_matrix<Matrix, SymmGroup> tmp;
                MPSTensor<Matrix, SymmGroup> & mps = mpsref.get();
                // prepare the LHS
                // Post-multiply Q^\dagger M on the std::vector<MPS<Matrix, SymmGroup> > elements
                // mps.make_left_paired(); // really?
                gemm(mps.data(), adjoint(Q), tmp);
                ret.push_back(tmp);

                // Now that we're ready with RHS, replace all MPSTensors with Q
                mps.data() = Q;
                mps.left_i = mps.data().left_basis();
                assert(mps.left_i == L.right_basis());
                mps.cur_normalization = Rnorm;
            }

            return ret;
        }
        else // if (method == SVD)
        {
            throw std::runtime_error("Multi_normalize_right with SVD NYI. Too lazy!!!");
        }

    }

    // assuming all elements in mps_vec have the same basis!!!
    return std::vector<block_matrix<Matrix, SymmGroup> >(mps_vec.size(),identity_matrix<block_matrix<Matrix, SymmGroup> >(mps_vec[0].get().data().left_basis()));
}

// Moves normalisation of an MPS vector to the right, just as the single-mps version in mps.hpp
template<class Matrix, class SymmGroup>
void multi_move_normalization_l2r(std::vector<MPS<Matrix, SymmGroup> > & vec, std::size_t p1, std::size_t p2, DecompMethod method)
{
    std::size_t length = vec[0].length();
    parallel::scheduler_balanced scheduler(length);

    // make sure all MPS are canonized at the same site
    std::size_t tmp_i = vec[0].canonized_i;

    for (auto&& mps: vec)
        assert(mps.canonized_i == tmp_i);
    for (int i = p1; i < std::min(p2, length); ++i)
    {
        Indicator normalization = vec[0][i].cur_normalization;
        for (auto&& mps: vec)
            assert(mps[i].cur_normalization == normalization);
        if (normalization == Lnorm) continue;

        // prepare a vector of MPSTensors for multi_normalize_left()
        std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > > mpst_vec;
        for (auto&& mps: vec)
            mpst_vec.push_back(std::ref(mps[i]));

        // right-hand side
        std::vector<block_matrix<Matrix, SymmGroup> > rhs;
        {
            parallel::guard proc(scheduler(i));
            rhs = multi_normalize_left(mpst_vec, method);
        }
        if (i < length-1) {
            parallel::guard proc(scheduler(i+1));

            assert(vec.size() == rhs.size());

            // TODO: transform this piece of code to fancy C++11
            // using std::transform and lambdas etc.

            for (int s = 0; s < vec.size(); s++)
            {
                vec[s][i+1].multiply_from_left(rhs[s]);
                vec[s][i+1].divide_by_scalar(vec[s][i+1].scalar_norm());
            }
        }
    }

    for (auto&& mps: vec)
        if (tmp_i == p1)
            mps.canonized_i = p2;
        else
            mps.canonized_i = std::numeric_limits<std::size_t>::max();
}

// Moves normalisation of an MPS vector to the right, just as the single-mps version in mps.hpp
template<class Matrix, class SymmGroup>
void multi_move_normalization_r2l(std::vector<MPS<Matrix, SymmGroup> > & vec, std::size_t p1, std::size_t p2, DecompMethod method)
{
    std::size_t length = vec[0].length();
    parallel::scheduler_balanced scheduler(length);

    // make sure all MPS are canonized at the same site
    std::size_t tmp_i = vec[0].canonized_i;

    for (auto&& mps: vec)
        assert(mps.canonized_i == tmp_i);

    for (int i = p1; i > std::max(p2, std::size_t(0)); --i)
    {
        Indicator normalization = vec[0][i].cur_normalization;
        for (auto&& mps: vec)
            assert(mps[i].cur_normalization == normalization);
        if (normalization == Rnorm) continue;

        // prepare a vector of MPSTensors for multi_normalize_left()
        std::vector<std::reference_wrapper<MPSTensor<Matrix, SymmGroup> > > mpst_vec;
        for (auto&& mps: vec)
            mpst_vec.push_back(std::ref(mps[i]));

        // left-hand side
        std::vector<block_matrix<Matrix, SymmGroup> > lhs;
        {
            parallel::guard proc(scheduler(i));
            lhs = multi_normalize_right(mpst_vec, method);
        }
        if (i > 0) {
            parallel::guard proc(scheduler(i-1));

            assert(vec.size() == lhs.size());

            // TODO: transform this piece of code to fancy new C++
            // using std::transform and lambdas etc.

            for (int s = 0; s < vec.size(); s++)
            {
                vec[s][i-1].multiply_from_right(lhs[s]);
                vec[s][i-1].divide_by_scalar(vec[s][i-1].scalar_norm());
            }
        }
    }

    for (auto&& mps: vec)
        if (tmp_i == p1)
            mps.canonized_i = p2;
        else
            mps.canonized_i = std::numeric_limits<std::size_t>::max();
}

template<class Matrix, class SymmGroup>
void multi_canonize(std::vector<MPS<Matrix, SymmGroup> > & vec, std::size_t center, DecompMethod method)
{
    std::size_t length = vec[0].length();
    std::size_t canonized_i = vec[0].canonized_i;
    for (auto&& mps: vec)
        assert(mps.canonized_i == canonized_i);

    if (canonized_i == center)
        return;

    if (canonized_i < center)
        multi_move_normalization_l2r(vec, canonized_i, center, method);
    else if (canonized_i < length)
        multi_move_normalization_r2l(vec, canonized_i, center, method);
    else {
        multi_move_normalization_l2r(vec, 0, center, method);
        multi_move_normalization_r2l(vec, length-1, center, method);
    }
    for (auto&& mps: vec)
        mps.canonized_i = center;
}

#endif