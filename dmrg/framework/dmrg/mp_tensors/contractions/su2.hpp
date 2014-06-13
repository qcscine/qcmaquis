#ifndef SU2_HPP
#define SU2_HPP

#include <alps/numeric/matrix.hpp>
#include <alps/numeric/matrix/algorithms.hpp>
#include "dmrg/block_matrix/detail/alps.hpp"
#include "dmrg/mp_tensors/mps.h"
#include "dmrg/mp_tensors/mps_mpo_ops.h"

#include <boost/tuple/tuple_io.hpp>

extern "C" {
    double gsl_sf_coupling_6j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf);
    double gsl_sf_coupling_9j(int two_ja, int two_jb, int two_jc, int two_jd, int two_je, int two_jf, int two_jg, int two_jh, int two_ji);
}

namespace SU2 {

    double mod_coupling(int two_ja, int two_jb, int two_jc,
                        int two_jd, int two_je, int two_jf,
                        int two_jg, int two_jh, int two_ji)
    {
        return sqrt( (two_jg+1.) * (two_jh+1.) * (two_jc+1.) * (two_jf+1.) )
        //return 1.0
                * gsl_sf_coupling_9j(two_ja, two_jb, two_jc,
                                     two_jd, two_je, two_jf,
                                     two_jg, two_jh, two_ji);
    }

    template<class Matrix1, class Matrix2, class Matrix3, class SymmGroup>
    void gemm(block_matrix<Matrix1, SymmGroup> const & A,
              block_matrix<Matrix2, SymmGroup> const & B,
              block_matrix<Matrix3, SymmGroup> & C)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Index<SymmGroup>::const_iterator ci;

        C.clear();
        assert(B.basis().is_sorted());

        const_iterator BBbegin = B.basis().begin();
        for (std::size_t k = 0; k < A.n_blocks(); ++k) {

            std::pair<const_iterator, const_iterator>
              er = std::equal_range(BBbegin, B.basis().end(),
                boost::make_tuple(A.right_basis_charge(k), SymmGroup::IdentityCharge, 0, 0), dual_index_detail::gt_row<SymmGroup>());

            for (const_iterator it = er.first; it != er.second; ++it)
            {
                std::size_t matched_block = std::distance(BBbegin, it);
                assert(matched_block == B.left_basis().position(A.right_basis_charge(k)));
                std::size_t new_block = C.insert_block(new Matrix3(num_rows(A[k]), num_cols(B[matched_block])),
                                                   A.left_basis_charge(k), B.right_basis_charge(matched_block));
                gemm(A[k], B[matched_block], C[new_block]);
            }
        }
    }
}

namespace contraction {
namespace SU2 {

    template<class Matrix, class OtherMatrix, class SymmGroup>
    void lbtm_kernel(size_t b2,
                     ContractionGrid<Matrix, SymmGroup>& contr_grid,
                     Boundary<OtherMatrix, SymmGroup> const & left,
                     std::vector<block_matrix<Matrix, SymmGroup> > const & left_mult_mps,
                     MPOTensor<Matrix, SymmGroup> const & mpo,
                     Index<SymmGroup> const & physical_i,
                     Index<SymmGroup> const & right_i,
                     Index<SymmGroup> const & out_left_i,
                     ProductBasis<SymmGroup> const & in_right_pb,
                     ProductBasis<SymmGroup> const & out_left_pb)
    {
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::index_type index_type;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::row_proxy row_proxy;
        typedef typename MPOTensor<OtherMatrix, SymmGroup>::col_proxy col_proxy;

        typedef typename SymmGroup::charge charge;
        typedef std::size_t size_t;

        col_proxy col_b2 = mpo.column(b2);
        for (typename col_proxy::const_iterator col_it = col_b2.begin(); col_it != col_b2.end(); ++col_it) {
            index_type b1 = col_it.index();

            block_matrix<Matrix, SymmGroup> const & T = left_mult_mps[b1];
            if (T.n_blocks() == 0) continue;
            MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(b1,b2);
            block_matrix<Matrix, SymmGroup> const & W = access.op;
            if (W.n_blocks() == 0) continue;

            // charge deltas are constant for all blocks
            charge operator_delta = SymmGroup::fuse(W.right_basis_charge(0), -W.left_basis_charge(0));
            charge        T_delta = SymmGroup::fuse(T.right_basis_charge(0), -T.left_basis_charge(0));
            charge    total_delta = SymmGroup::fuse(operator_delta, -T_delta);

            block_matrix<Matrix, SymmGroup>& ret = contr_grid(b1,b2);

            for (size_t r = 0; r < right_i.size(); ++r)
            {
                charge out_r_charge = right_i[r].first;
                charge out_l_charge = SymmGroup::fuse(out_r_charge, total_delta);
                size_t r_size = right_i[r].second;

                if (!out_left_i.has(out_l_charge)) continue;

                size_t o = ret.find_block(out_l_charge, out_r_charge);
                if ( o == ret.n_blocks() ) {
                    o = ret.insert_block(Matrix(1,1), out_l_charge, out_r_charge);
                    ret.resize_block(out_l_charge, out_r_charge, out_left_i.size_of_block(out_l_charge), r_size);
                }

                for (size_t w_block = 0; w_block < W.n_blocks(); ++w_block)
                {
                    charge phys_c1 = W.left_basis_charge(w_block);
                    charge phys_c2 = W.right_basis_charge(w_block);

                    charge in_r_charge = SymmGroup::fuse(out_r_charge, -phys_c1);
                    charge in_l_charge = SymmGroup::fuse(in_r_charge, -T_delta);
                    size_t t_block = T.basis().position(in_l_charge, in_r_charge);
                    if (t_block == T.n_blocks()) continue;

                    size_t in_right_offset = in_right_pb(phys_c1, out_r_charge);
                    size_t out_left_offset = out_left_pb(phys_c2, in_l_charge);

                    size_t phys_s1 = W.left_basis_size(w_block);
                    size_t phys_s2 = W.right_basis_size(w_block);
                    Matrix const & wblock = W[w_block];
                    Matrix const & iblock = T[t_block];
                    Matrix & oblock = ret[o];

                    maquis::dmrg::detail::lb_tensor_mpo(oblock, iblock, wblock,
                            out_left_offset, in_right_offset,
                            phys_s1, phys_s2, T.left_basis_size(t_block), r_size, access.scale);
                }
            } // right index block
        } // b1
    }

    template <class SymmGroup>
    bool column_check(typename SymmGroup::charge c, Index<SymmGroup> const & ind)
    {
        int count = 0;
        for (std::size_t p = 0; p < ind.size(); ++p)
            if (ind[p].first == c) ++count;

        return (count == 1);
    }

    double couple_destroy(int jR, int jRt, int local_spin)
    {
        // jR = bra right spin, jRt = ket right spin

        double ret = 1.;
        if (local_spin == 0) {
            ret = std::sqrt( (jR + 1.) / (jRt + 1.));
            int phase = (jRt - jR + 1)/2;
            for (int p_= 0; p_ < phase; ++p_) ret *= -1.;
        }
        return ret;
    }

    double couple_create(int jR, int jLt, int local_spin, int spin_out)
    {
        // jR = bra..
        double ret;
        if(local_spin == 0) {
            ret = 1./sqrt(2.);
        }
        else {
            ret = std::sqrt( (jLt + 1.) / (jR + 1.));
            ret *= 1./sqrt(2.);
            int phase = (jR - jLt + 1)/2;
            for (int p_= 0; p_ < phase; ++p_) ret *= -1.;
        }
        return ret;
    }

    template <typename T> int sgn(T val) {
        return (T(0) < val) - (val < T(0));
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    block_matrix<OtherMatrix, SymmGroup>
    apply_operator(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          block_matrix<OtherMatrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo,
                          int boundary_spin, std::vector<int> const & config,
                          bool debug = false)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;
        typedef typename DualIndex<SymmGroup>::const_iterator const_iterator;
        typedef typename Index<SymmGroup>::const_iterator ci;

        assert(ket_tensor.phys_i == bra_tensor.phys_i);

        bra_tensor.make_left_paired();
        ket_tensor.make_right_paired();

        MPOTensor_detail::const_term_descriptor<Matrix, SymmGroup> access = mpo.at(0,0);
        block_matrix<Matrix, SymmGroup> const & W = access.op;
        //if(debug) maquis::cout << W << std::endl;
        block_matrix<OtherMatrix, SymmGroup> ret;

        block_matrix<OtherMatrix, SymmGroup> t1;
        ::SU2::gemm(left, ket_tensor.data(), t1);


        Index<SymmGroup> const & left_i = ket_tensor.row_dim();
        Index<SymmGroup> const & right_i = ket_tensor.col_dim();
        Index<SymmGroup> const & phys_i = ket_tensor.site_dim();

        ProductBasis<SymmGroup> out_left_pb(phys_i, left_i);    
        ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));
        charge opdelta[2];
        opdelta[0] = phys_i[1].first; opdelta[1] = phys_i[2].first;

        for (std::size_t k = 0; k < left.n_blocks(); ++k) {

            std::pair<const_iterator, const_iterator>
              er = std::equal_range(ket_tensor.data().basis().begin(), ket_tensor.data().basis().end(),
                boost::make_tuple(left.right_basis_charge(k), SymmGroup::IdentityCharge, 0, 0),
                  dual_index_detail::gt_row<SymmGroup>());

            for (const_iterator it = er.first; it != er.second; ++it)
            {
                std::size_t matched_block = std::distance(ket_tensor.data().basis().begin(), it);
                //Matrix T_block(num_rows(left[k]), num_cols(ket_tensor[matched_block]));
                //gemm(A[k], B[matched_block], T_block);
                charge lc = left.left_basis_charge(k);
                charge mc = left.right_basis_charge(k);
                charge rc = ket_tensor.data().right_basis_charge(matched_block);
                charge mc1 = ket_tensor.data().left_basis_charge(matched_block);
                assert (mc == mc1);
                size_t t_pos = t1.basis().position(lc, rc);

                for (size_t w_block = 0; w_block < phys_i.size(); ++w_block)
                {
                    charge phys_in = W.left_basis_charge(w_block);
                    charge phys_out = W.right_basis_charge(w_block);

                    charge free_rc = SymmGroup::fuse(rc, phys_in);
                    if (!right_i.has(free_rc))
                        continue;

                    charge new_rc = SymmGroup::fuse(lc, phys_out);
                    if (!bra_tensor.col_dim().has(new_rc))
                        continue;


                    int i  = lc[1], ip = new_rc[1];
                    int j  = mc[1], jp  = free_rc[1];
                    int two_sp = std::abs(i - ip), two_s  = std::abs(j - jp);
                    int a = boundary_spin, k = std::abs(phys_in[1]-phys_out[1]), ap = (a + k == 2) ? 0 : a + k;
                    double coupling_coeff = ::SU2::mod_coupling(j, two_s, jp, a,k,ap, i, two_sp, ip);
                    coupling_coeff *= access.scale * W[w_block](0,0);

                    // Coupling coefficient
                    double coupling = access.scale * W[w_block](0,0);
                    if (SymmGroup::fuse(phys_out, -phys_in)[0] == -1) { // if destructor
                        coupling *= couple_destroy(new_rc[1], free_rc[1], phys_in[1]);
                    }
                    else if (SymmGroup::fuse(phys_out, -phys_in)[0] == 1) { // if creator
                        coupling *= couple_create(free_rc[1], mc[1], phys_in[1], 0);

                        // avoid Spin 1 Tensor
                        if (std::abs(new_rc[1]-free_rc[1]) == 2)
                            continue;
                    }
                    else if (std::abs(new_rc[1]-free_rc[1]) == 1) { // if Spin 1/2 tensor
                        assert(phys_in==phys_out);
                        if (std::abs(free_rc[0] - mc[0] == 1))
                        {
                            int phase = ((((j+ip)/2)%2)!=0)?-1:1;
                            coupling *= phase * sqrt((j+1.) * (ip+1.)) * gsl_sf_coupling_6j(ip, jp, 1, j, i, 1);
                        }

                        //if(config[4]) coupling *= ((((i - j + 3)/2)%2)!=0)?-1:1;
                        //if(config[5]) coupling *= ((((ip - jp + 3)/2)%2)!=0)?-1:1;
                        //coupling *= pow(ip+1., double(config[0])/2.) * pow(jp+1., double(config[1])/2.);
                        //coupling *= pow(i+1., double(config[2])/2.) * pow(j+1., double(config[3])/2.);
                        //coupling *= pow(two_s+1., double(config[6])/2.); //* pow(two_sp+1., double(config[7])/2.);
                        //coupling *= pow(a, double(config[8])/2.) * pow(ap, double(config[9])/2.);
                        //coupling *= pow(k, double(config[10])/2.)
                        //coupling *= coupling_coeff;
                    }
                    //coupling_coeff *= sgn(coupling_coeff) * sgn(coupling);

                    if (debug && std::abs(new_rc[1]-free_rc[1]) == 1) {
                        std::cout << j << "," << two_s << "," << jp << " | " << a << "," << k << "," << ap << " | "
                                  << i << "," << two_sp << "," << ip << " | " << phys_in << phys_out
                                  << std::right << std::setw(8) << "cc: " << std::setw(12) << coupling_coeff
                                  << " | " << std::setw(12) << coupling << std::endl;
                    }

                    // T Access
                    if (debug && std::abs(new_rc[1]-free_rc[1]) && a==1)
                    maquis::cout << "access " << mc << " + " << phys_in<< "|" << free_rc << " -- "
                                 << lc << " + " << phys_out << "|" << new_rc << std::endl;
                    size_t right_offset = in_right_pb(phys_in, free_rc);
                    //size_t ldim = left_i.size_of_block(mc);
                    size_t ldim = t1.left_basis_size(t_pos);
                    size_t rdim = right_i.size_of_block(free_rc);
                    //if(debug) maquis::cout << "t1 block@" << t_pos << " " << ldim << " " << rdim  << " " << right_offset << std::endl;
                    Matrix T_cp(ldim, rdim, 0);
                    for (size_t c=0; c<rdim; ++c)
                        std::transform(t1[t_pos].col(right_offset+c).first,
                                       t1[t_pos].col(right_offset+c).second, T_cp.col(c).first, boost::lambda::_1*coupling);

                    // Bra Access
                    size_t bra_index = bra_tensor.col_dim().position(new_rc);
                    assert(column_check(new_rc, bra_tensor.col_dim()));
                    Matrix const & bra_block = bra_tensor.data()[bra_index];
                    // mc + phys_out = new_rc
                    size_t left_offset = out_left_pb(phys_out, lc);
                    ldim = bra_tensor.row_dim().size_of_block(lc);
                    rdim = bra_tensor.col_dim()[bra_index].second;
                    //if(debug) maquis::cout << "brablock@" << bra_index << " " << ldim << " " << rdim  << " " << left_offset << std::endl;
                    Matrix bra_cp(ldim, rdim, 0);
                    for (size_t row=0; row<ldim; ++row)
                        std::copy(bra_block.row(left_offset+row).first,
                                  bra_block.row(left_offset+row).second, bra_cp.row(row).first);

                    // Multiply
                    Matrix prod(num_cols(bra_cp), num_cols(T_cp), 0);
                    gemm(transpose(bra_cp), T_cp, prod);

                    // Check-in
                    ret.match_and_add_block(prod, new_rc, free_rc);
                }
            }
        }
        if(debug) maquis::cout << std::endl;
        return ret;
    }

    template<class Matrix, class OtherMatrix, class SymmGroup>
    Boundary<OtherMatrix, SymmGroup>
    overlap_mpo_left_step(MPSTensor<Matrix, SymmGroup> const & bra_tensor,
                          MPSTensor<Matrix, SymmGroup> const & ket_tensor,
                          Boundary<OtherMatrix, SymmGroup> const & left,
                          MPOTensor<Matrix, SymmGroup> const & mpo)
    {
        typedef typename SymmGroup::charge charge;
        typedef typename MPOTensor<Matrix, SymmGroup>::index_type index_type;

        MPSTensor<Matrix, SymmGroup> ket_cpy = ket_tensor;
        std::vector<block_matrix<Matrix, SymmGroup> > t = boundary_times_mps(ket_cpy, left, mpo);

        Index<SymmGroup> const & left_i = bra_tensor.row_dim();
        Index<SymmGroup> const & right_i = ket_tensor.col_dim();
        Index<SymmGroup> out_left_i = ket_tensor.site_dim() * left_i;
        ProductBasis<SymmGroup> out_left_pb(ket_tensor.site_dim(), left_i);
        ProductBasis<SymmGroup> in_right_pb(ket_tensor.site_dim(), right_i,
                                boost::lambda::bind(static_cast<charge(*)(charge, charge)>(SymmGroup::fuse),
                                        -boost::lambda::_1, boost::lambda::_2));

        index_type loop_max = mpo.col_dim();

        bra_tensor.make_left_paired();
        block_matrix<Matrix, SymmGroup> bra_conj = conjugate(bra_tensor.data());

        Boundary<Matrix, SymmGroup> ret;
        ret.resize(loop_max);

        omp_for(index_type b2, range<index_type>(0,loop_max), {
            ContractionGrid<Matrix, SymmGroup> contr_grid(mpo, 0, 0);
            //maquis::cout << t[b2] << std::endl;
            contraction::SU2::lbtm_kernel(b2, contr_grid, left, t, mpo, ket_tensor.site_dim(), right_i, out_left_i, in_right_pb, out_left_pb);
            //maquis::cout << contr_grid(0,0) << "---------------------" << std::endl << std::endl;
            ::SU2::gemm(transpose(contr_grid(0,0)), bra_conj, ret[b2]);
        });

        return ret;
    }

} // namespace SU2
} // namespace contraction

namespace SU2 {

    template<class Matrix, class SymmGroup>
    double expval(MPS<Matrix, SymmGroup> const & mps, MPO<Matrix, SymmGroup> const & mpo,
                  int p1, int p2, std::vector<int> config)
    {
        bool debug = false;
        if (p1 == 1 && p2 == 3) debug = true;

        assert(mpo.length() == mps.length());
        std::size_t L = mps.length();
        Boundary<Matrix, SymmGroup> left = mps.left_boundary();

        for(size_t i = 0; i < L; ++i) {
            MPSTensor<Matrix, SymmGroup> cpy = mps[i];
            if (i==p1) 
                left[0] = contraction::SU2::apply_operator(cpy, mps[i], left[0], mpo[i], 0, config, debug);
            else if (p1 < i && i < p2)
                left[0] = contraction::SU2::apply_operator(cpy, mps[i], left[0], mpo[i], 1, config, debug);
            else if (i==p2) 
                left[0] = contraction::SU2::apply_operator(cpy, mps[i], left[0], mpo[i], 1, config, debug);
            else 
                left[0] = contraction::SU2::apply_operator(cpy, mps[i], left[0], mpo[i], 0, config, debug);
                //left = contraction::SU2::overlap_mpo_left_step(mps[i], mps[i], left, mpo[i]);

            //if (i==1 && p1 == 1 && p2 == 2) maquis::cout << left[0] << std::endl;
            //if (i==2 && p1 == 1 && p2 == 2) maquis::cout << left[0] << std::endl;
        }

        return maquis::real(left[0].trace());
    }

} // namespace SU2

#endif
