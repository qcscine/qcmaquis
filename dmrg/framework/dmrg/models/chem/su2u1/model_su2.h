/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
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

#ifndef QC_SU2_H
#define QC_SU2_H

#include <cmath>
#include <sstream>
#include <fstream>
#include <iterator>
#include <boost/shared_ptr.hpp>
#include <boost/tokenizer.hpp>
#include <boost/regex.hpp>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

#include "dmrg/models/chem/util.h"
#include "dmrg/models/chem/parse_integrals.h"
#include "dmrg/models/chem/pg_util.h"
#include "dmrg/models/chem/su2u1/chem_helper_su2.h"
#include "dmrg/models/chem/su2u1/operator_spin_variants.h"
#include "dmrg/models/chem/su2u1/term_maker_su2.h"

template<class Matrix, class SymmGroup>
class qc_su2 : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Lattice::pos_t pos_t;
    typedef typename Matrix::value_type value_type;
    typedef typename alps::numeric::associated_one_matrix<Matrix>::type one_matrix;

public:
    
    qc_su2(Lattice const & lat_, BaseParameters & parms_);
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    // For this model: site_type == point group irrep
    Index<SymmGroup> const & phys_dim(size_t type) const
    {
        return phys_indices[type];
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill;
    }

    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms_) const
    {
        return chem_detail::qn_helper<SymmGroup>().total_qn(parms_);
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "docc")
            return docc;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    measurements_type measurements () const
    {
        measurements_type meas;
        return meas;
    }

private:
    Lattice const & lat;
    BaseParameters & parms;
    std::vector<Index<SymmGroup> > phys_indices;

    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    tag_type create_fill, create, destroy_fill, destroy,
             create_fill_count, create_count, destroy_fill_count, destroy_count,
             count, docc, e2d, d2e, flip_S0, flip_to_S2, flip_to_S0,
             ident, ident_full, fill, count_fill;

    typename SymmGroup::subcharge max_irrep;

    std::vector<op_t> generate_site_specific_ops(op_t const & op) const
    {
        PGDecorator<SymmGroup> set_symm;
        std::vector<op_t> ret;
        for (typename SymmGroup::subcharge sc=0; sc < max_irrep+1; ++sc) {
            op_t mod(set_symm(op.left_basis(), sc), set_symm(op.right_basis(), sc));
            for (std::size_t b = 0; b < op.n_blocks(); ++b)
                mod[b] = op[b];

            ret.push_back(mod);
        } return ret;
    }

};

template <class Matrix, class SymmGroup>
qc_su2<Matrix, SymmGroup>::qc_su2(Lattice const & lat_, BaseParameters & parms_)
: lat(lat_)
, parms(parms_)
, tag_handler(new table_type())
{
    typedef typename SymmGroup::subcharge subcharge;
    // find the highest irreducible representation number 
    // used to generate ops for all irreps 0..max_irrep
    max_irrep = 0;
    for (pos_t p=0; p < lat.size(); ++p)
        max_irrep = (lat.get_prop<typename SymmGroup::subcharge>("type", p) > max_irrep)
                        ? lat.get_prop<typename SymmGroup::subcharge>("type", p) : max_irrep;

    typename SymmGroup::charge A(0), B(0), C(0), D(0);
    A[0] = 2; // 20
    B[0] = 1; B[1] =  1; // 11
    C[0] = 1; C[1] = -1; // 1-1
    // D = 00

    SpinDescriptor<symm_traits::SU2Tag> one_half_up(1,0,1);
    SpinDescriptor<symm_traits::SU2Tag> one_half_down(1,1,0);
    SpinDescriptor<symm_traits::SU2Tag> one_up(2,0,2);
    SpinDescriptor<symm_traits::SU2Tag> one_flat(2,1,1);
    SpinDescriptor<symm_traits::SU2Tag> one_down(2,2,0);

    for (subcharge irr=0; irr <= max_irrep; ++irr)
    {
        Index<SymmGroup> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(B, irr), 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(C, irr), 1));
        phys.insert(std::make_pair(D, 1));

        phys_indices.push_back(phys);
    }

    // cheaper to use this for spin0 tensors, instead of ident_full
    op_t ident_op;
    ident_op.insert_block(Matrix(1,1,1), A, A);
    ident_op.insert_block(Matrix(1,1,1), B, B);
    ident_op.insert_block(Matrix(1,1,1), C, C);
    ident_op.insert_block(Matrix(1,1,1), D, D);

    // apply if spin > 0
    op_t ident_full_op;
    ident_full_op.insert_block(Matrix(1,1,1), A, A);
    ident_full_op.insert_block(Matrix(1,1,1), D, D);
    ident_full_op.insert_block(Matrix(1,1,1), B, B);
    ident_full_op.insert_block(Matrix(1,1,1), C, C);
    ident_full_op.insert_block(Matrix(1,1,1), B, C);
    ident_full_op.insert_block(Matrix(1,1,1), C, B);

    op_t fill_op;
    fill_op.insert_block(Matrix(1,1,1),  A, A);
    fill_op.insert_block(Matrix(1,1,1),  D, D);
    fill_op.insert_block(Matrix(1,1,-1), B, B);
    fill_op.insert_block(Matrix(1,1,-1), C, C);
    fill_op.insert_block(Matrix(1,1,-1), B, C);
    fill_op.insert_block(Matrix(1,1,-1), C, B);

    /*************************************************************/

    op_t create_fill_op;
    create_fill_op.spin = one_half_up;
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_fill_op.insert_block(Matrix(1,1,1), D, B);
    create_fill_op.insert_block(Matrix(1,1,1), D, C);

    op_t destroy_op;
    destroy_op.spin = one_half_down;
    destroy_op.insert_block(Matrix(1,1,1), A, B);
    destroy_op.insert_block(Matrix(1,1,1), A, C);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

    op_t destroy_fill_op;
    destroy_fill_op.spin = one_half_up;
    destroy_fill_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_op.insert_block(Matrix(1,1,1), A, C);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), B, D);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

    op_t create_op;
    create_op.spin = one_half_down;
    create_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_op.insert_block(Matrix(1,1,-1), D, B);
    create_op.insert_block(Matrix(1,1,-1), D, C);

    /*************************************************************/

    op_t create_fill_count_op;
    create_fill_count_op.spin = one_half_up;
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

    op_t destroy_count_op;
    destroy_count_op.spin = one_half_down;
    destroy_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t destroy_fill_count_op;
    destroy_fill_count_op.spin = one_half_up;
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t create_count_op;
    create_count_op.spin = one_half_down;
    create_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

    /*************************************************************/

    op_t count_op;
    count_op.insert_block(Matrix(1,1,2), A, A);
    count_op.insert_block(Matrix(1,1,1), B, B);
    count_op.insert_block(Matrix(1,1,1), C, C);

    op_t docc_op;
    docc_op.insert_block(Matrix(1,1,1), A, A);

    op_t e2d_op;
    e2d_op.insert_block(Matrix(1,1,1), D, A);

    op_t d2e_op;
    d2e_op.insert_block(Matrix(1,1,1), A, D);

    op_t count_fill_op;
    count_fill_op.insert_block(Matrix(1,1,2),  A, A);
    count_fill_op.insert_block(Matrix(1,1,-1), B, B);
    count_fill_op.insert_block(Matrix(1,1,-1), C, C);
    count_fill_op.insert_block(Matrix(1,1,-1), B, C);
    count_fill_op.insert_block(Matrix(1,1,-1), C, B);

    op_t flip_to_S2_op;
    flip_to_S2_op.spin = one_up;
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2)), B, B);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)), C, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  B, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  C, B);

    op_t flip_to_S0_op = flip_to_S2_op;
    flip_to_S0_op.spin = one_down;

    op_t flip_S0_op = flip_to_S2_op;
    flip_S0_op.spin = one_flat;

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);

    REGISTER(ident,        tag_detail::bosonic)
    REGISTER(ident_full,   tag_detail::bosonic)
    REGISTER(fill,         tag_detail::bosonic)

    REGISTER(create_fill,  tag_detail::fermionic)
    REGISTER(create,       tag_detail::fermionic)
    REGISTER(destroy_fill, tag_detail::fermionic)
    REGISTER(destroy,      tag_detail::fermionic)

    REGISTER(create_fill_count,  tag_detail::fermionic)
    REGISTER(create_count,       tag_detail::fermionic)
    REGISTER(destroy_fill_count, tag_detail::fermionic)
    REGISTER(destroy_count,      tag_detail::fermionic)

    REGISTER(count,        tag_detail::bosonic)
    REGISTER(docc,         tag_detail::bosonic)
    REGISTER(e2d,          tag_detail::bosonic)
    REGISTER(d2e,          tag_detail::bosonic)
    REGISTER(flip_S0,      tag_detail::bosonic)
    REGISTER(flip_to_S2,   tag_detail::bosonic)
    REGISTER(flip_to_S0,   tag_detail::bosonic)
    REGISTER(count_fill,   tag_detail::bosonic)

#undef REGISTER

//#define PRINT(op) maquis::cout << #op << "\t" << op << std::endl;
//    PRINT(ident)
//    PRINT(ident_full)
//    PRINT(fill)
//    PRINT(create_fill)
//    PRINT(create)
//    PRINT(destroy_fill)
//    PRINT(destroy)
//    PRINT(count)
//    PRINT(count_fill)
//    PRINT(docc)
//    PRINT(e2d)
//    PRINT(d2e)
//#undef PRINT

    /**********************************************************************/
    OperatorSpinVariants<Matrix, SymmGroup> creators(create, tag_handler);
    OperatorSpinVariants<Matrix, SymmGroup> creators_fill(create_fill, tag_handler);
    OperatorSpinVariants<Matrix, SymmGroup> destructors(destroy, tag_handler);
    OperatorSpinVariants<Matrix, SymmGroup> destructors_fill(destroy_fill, tag_handler);

    typename TermMakerSU2<Matrix, SymmGroup>::OperatorPackage create_pkg(&creators, &creators_fill)
                                                            , destroy_pkg(&destructors, &destructors_fill);
    /**********************************************************************/

    chem_detail::ChemHelperSU2<Matrix, SymmGroup> ta(parms, lat, ident, ident, tag_handler);
    alps::numeric::matrix<Lattice::pos_t> idx_ = ta.getIdx();
    std::vector<value_type> matrix_elements = ta.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);
 
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = idx_(m, 0);
        int j = idx_(m, 1);
        int k = idx_(m, 2);
        int l = idx_(m, 3);

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back( boost::make_tuple(0, ident) );
            this->terms_.push_back(term);
            
            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back( boost::make_tuple(i, count));
            this->terms_.push_back(term);

            used_elements[m] += 1;
            continue;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {

            // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, ident, std::sqrt(2.)*matrix_elements[m],i,j,create, create_fill, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, ident, std::sqrt(2.)*matrix_elements[m],j,i,create, create_fill, destroy, destroy_fill
            ));

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back(boost::make_tuple(i, docc));
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {

            int same_idx, pos1;

            if      (i==j) { same_idx = i; pos1 = l; }
            else if (k==l) { same_idx = l; pos1 = i; }
            else           { throw std::runtime_error("Term generation logic has failed for V_ijjj term\n"); }

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, ident,  std::sqrt(2.)*matrix_elements[m], same_idx, pos1, create_count, create_fill_count, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, ident, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, destroy_count, destroy_fill_count, create, create_fill
            ));

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            ta.add_2term(this->terms_, TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, k, count, count));

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, j, e2d, d2e));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, j, d2e, e2d));

            // here we have spin0--j--spin1--i--spin0
            // the sqrt(3.) counteracts the Clebsch coeff C^{110}_{mrm'} which applies when the spin1 couples back to spin0
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                false, ident_full, std::sqrt(3.) * matrix_elements[m], i, j, flip_to_S0, flip_to_S2, flip_to_S0, flip_to_S2
            ));

            ta.add_2term(this->terms_, TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, -0.5 * matrix_elements[m], i, j, count, count));

            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {
            typedef TermMakerSU2<Matrix, SymmGroup> TM;

            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; k = i; l = j; }

            std::vector<term_descriptor> & vec = this->terms_;

            ta.add_3term(vec, TM::three_term(ident, std::sqrt(2.)*matrix_elements[m], same_idx, k, l,
                                             count, count_fill, create, create_fill, destroy, destroy_fill));
            ta.add_3term(vec, TM::three_term(ident, std::sqrt(2.)*matrix_elements[m], same_idx, l, k,
                                             count, count_fill, create, create_fill, destroy, destroy_fill));

            used_elements[m] += 1;
        }

        // 9887 7371 8727

        // 4-fold degenerate (+spin) V_ijil = V_ijli = V_jiil = V_jili  <--- coded
        //                           V_ilij = V_ilji = V_liij = V_liji
        else if ( ((i==k && j!=l) || j==k || (j==l && i!=k)) && (i!=j && k!=l)) {
            typedef TermMakerSU2<Matrix, SymmGroup> TM;

            int same_idx, pos1, pos2;
            if (i==k) { same_idx = i; pos1 = l; pos2 = j; }
            if (j==k) { same_idx = j; pos1 = l; pos2 = i; }
            if (j==l) { same_idx = j; pos1 = k; pos2 = i; }

            std::vector<term_descriptor> & vec = this->terms_;

            // Note: need minus because of clebsch gordan coeff from two destructors or two creators
            vec.push_back(TM::three_term(ident, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, e2d, e2d, destroy, destroy_fill, destroy, destroy_fill)); vec.push_back(TM::three_term(ident, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, d2e, d2e, create, create_fill, create, create_fill));

            if ( same_idx < std::min(pos1,pos2) )
            {
                this->terms_.push_back(TM::three_term(
                    ident_full, std::sqrt(3.)*matrix_elements[m], same_idx, pos1, pos2, flip_to_S2, flip_to_S2, create, creators_fill(2,1), destroy, destructors_fill(2,1)
                ));
                ta.add_3term(vec, TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, count, count, create, create_fill, destroy, destroy_fill
                ));

                this->terms_.push_back(TM::three_term(
                    // note minus sign, because commutation on same_idx is not taken into account
                    ident_full, -std::sqrt(3.)*matrix_elements[m], same_idx, pos2, pos1, flip_to_S2, flip_to_S2, create, creators_fill(2,1), destroy, destructors_fill(2,1)
                ));
                ta.add_3term(vec, TM::three_term(
                    ident,  -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos2, pos1, count, count, create, create_fill, destroy, destroy_fill
                ));
            }
            else if (same_idx > std::max(pos1,pos2))
            {
                this->terms_.push_back(TM::three_term(
                    ident_full, std::sqrt(3.)*matrix_elements[m], same_idx, pos1, pos2, flip_to_S0, flip_to_S0, creators(1,2), create_fill, destructors(1,2), destroy_fill
                ));
                ta.add_3term(vec, TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, count, count, create, create_fill, destroy, destroy_fill
                ));

                this->terms_.push_back(TM::three_term(
                    // note minus sign, because commutation on same_idx is not taken into account
                    ident_full, -std::sqrt(3.)*matrix_elements[m], same_idx, pos2, pos1, flip_to_S0, flip_to_S0, creators(1,2), create_fill, destructors(1,2), destroy_fill
                ));
                ta.add_3term(vec, TM::three_term(
                    ident,  -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos2, pos1, count, count, create, create_fill, destroy, destroy_fill
                ));
            }
            else
            {
                this->terms_.push_back(TM::three_term(
                    ident,      std::sqrt(3.)*matrix_elements[m], same_idx, pos1, pos2, flip_S0, flip_S0, create, create_fill, destroy, destroy_fill
                ));
                ta.add_3term(vec, TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, count_fill, count_fill, create, create_fill, destroy, destroy_fill
                ));

                this->terms_.push_back(TM::three_term(
                    ident,     -std::sqrt(3.)*matrix_elements[m], same_idx, pos2, pos1, flip_S0, flip_S0, create, create_fill, destroy, destroy_fill
                ));
                ta.add_3term(vec, TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos2, pos1, count_fill, count_fill, create, create_fill, destroy, destroy_fill
                ));
            }

            used_elements[m] += 1;
        }

        // 32 (8x4)-fold degenerate V_ijkl = V_jikl = V_ijlk = V_jilk = V_klij = V_lkij = V_klji = V_lkji * spin
        // V_ijkl -> 24 permutations which fall into 3 equivalence classes of 8 elements (with identical V_ matrix element)
        // coded: 4 index permutations which generate all Sz spins
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {
            typedef TermMakerSU2<Matrix, SymmGroup> TM;
            std::vector<term_descriptor> & vec = this->terms_;

            // These 3 cases produce different S_z spin patterns, which differ along with different index permutations
            // As in standard notation of the Hamiltonian, the first two positions get a creator, the last two a destructor

            if (k > l && l > j) // eg V_4132
            { // generates up|up|up|up, up|down|down|up, down|up|up|down, down|down|down|down
                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], i,k,l,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], i,k,l,j, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], i,l,k,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], i,l,k,j, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], j,k,l,i, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], j,k,l,i, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], j,l,k,i, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], j,l,k,i, create_pkg, destroy_pkg));
            }
            else if (k > j && j > l) // eg V_4231
            { // generates up|up|up|up, up|down|up|down, down|up|down|up, down|down|down|down
                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], i,k,l,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,               -matrix_elements[m], i,k,l,j, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2,  std::sqrt(3.)*matrix_elements[m], i,l,k,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], i,l,k,j, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2,  std::sqrt(3.)*matrix_elements[m], j,k,l,i, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,                matrix_elements[m], j,k,l,i, create_pkg, destroy_pkg));

                ta.add_4term(vec, TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_elements[m], j,l,k,i, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident,      1,               -matrix_elements[m], j,l,k,i, create_pkg, destroy_pkg));
            }
            else if (j > k && k > l) // eg V_4321
            { // generates up|up|up|up, up|up|down|down, down|down|up|up, down|down|down|down
                ta.add_4term(vec, TM::four_term(ident, 1, 2.*matrix_elements[m], i,k,l,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident, 1, 2.*matrix_elements[m], i,l,k,j, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident, 1, 2.*matrix_elements[m], j,k,l,i, create_pkg, destroy_pkg));
                ta.add_4term(vec, TM::four_term(ident, 1, 2.*matrix_elements[m], j,l,k,i, create_pkg, destroy_pkg));
            }
            else { throw std::runtime_error("unexpected index arrangment in V_ijkl term\n"); }

            used_elements[m] += 1;
        }

    } // matrix_elements for

    ta.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
