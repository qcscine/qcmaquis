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
             create_fill_couple_down, destroy_fill_couple_down,
             create_couple_up, destroy_couple_up,
             create_fill_count, create_count, destroy_fill_count, destroy_count,
             count, docc, e2d, d2e, flip_to_S2, flip_to_S0,
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
    create_fill_op.twoS = 1; create_fill_op.twoSaction = 1;
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_fill_op.insert_block(Matrix(1,1,1), D, B);
    create_fill_op.insert_block(Matrix(1,1,1), D, C);

    op_t destroy_op;
    destroy_op.twoS = 1; destroy_op.twoSaction = -1;
    destroy_op.insert_block(Matrix(1,1,1), A, B);
    destroy_op.insert_block(Matrix(1,1,1), A, C);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

    op_t destroy_fill_op;
    destroy_fill_op.twoS = 1; destroy_fill_op.twoSaction = 1;
    destroy_fill_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_op.insert_block(Matrix(1,1,1), A, C);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), B, D);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

    op_t create_op;
    create_op.twoS = 1; create_op.twoSaction = -1;
    create_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_op.insert_block(Matrix(1,1,-1), D, B);
    create_op.insert_block(Matrix(1,1,-1), D, C);

    /*************************************************************/

    op_t create_fill_couple_down_op = create_fill_op;
    create_fill_couple_down_op.twoSaction = -1;

    op_t destroy_fill_couple_down_op = destroy_fill_op;
    destroy_fill_couple_down_op.twoSaction = -1;

    op_t create_couple_up_op = create_op;
    create_couple_up_op.twoSaction = 1;

    op_t destroy_couple_up_op = destroy_op;
    destroy_couple_up_op.twoSaction = 1;

    /*************************************************************/

    op_t create_fill_count_op;
    create_fill_count_op.twoS = 1; create_fill_count_op.twoSaction = 1;
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

    op_t destroy_count_op;
    destroy_count_op.twoS = 1; destroy_count_op.twoSaction = -1;
    destroy_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t destroy_fill_count_op;
    destroy_fill_count_op.twoS = 1; destroy_fill_count_op.twoSaction = 1;
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t create_count_op;
    create_count_op.twoS = 1; create_count_op.twoSaction = -1;
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
    flip_to_S2_op.twoS = 2; flip_to_S2_op.twoSaction = 2;
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2)), B, B);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)), C, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  B, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  C, B);

    op_t flip_to_S0_op = flip_to_S2_op;
    flip_to_S0_op.twoSaction = -2;

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

    REGISTER(create_fill_couple_down,  tag_detail::fermionic)
    REGISTER(destroy_fill_couple_down,  tag_detail::fermionic)
    REGISTER(create_couple_up,  tag_detail::fermionic)
    REGISTER(destroy_couple_up,  tag_detail::fermionic)

    REGISTER(create_fill_count,  tag_detail::fermionic)
    REGISTER(create_count,       tag_detail::fermionic)
    REGISTER(destroy_fill_count, tag_detail::fermionic)
    REGISTER(destroy_count,      tag_detail::fermionic)

    REGISTER(count,        tag_detail::bosonic)
    REGISTER(docc,         tag_detail::bosonic)
    REGISTER(e2d,          tag_detail::bosonic)
    REGISTER(d2e,          tag_detail::bosonic)
    REGISTER(flip_to_S2,   tag_detail::bosonic)
    REGISTER(flip_to_S0,   tag_detail::bosonic)
    REGISTER(count_fill,   tag_detail::bosonic)

#undef REGISTER
    /**********************************************************************/

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

    //chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, ident, ident, tag_handler);
    alps::numeric::matrix<Lattice::pos_t> idx_;
    std::vector<value_type> matrix_elements;
    boost::tie(idx_, matrix_elements) = chem_detail::parse_integrals<value_type>(parms, lat);

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
            {
                term_descriptor term;
                term.coeff = matrix_elements[m];
                term.push_back( boost::make_tuple(i, count));
                this->terms_.push_back(term);
            }

            used_elements[m] += 1;
            continue;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {

            // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, lat.size(), ident, fill, std::sqrt(2.)*matrix_elements[m],i,j,create, create_fill, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, lat.size(), ident, fill, std::sqrt(2.)*matrix_elements[m],j,i,create, create_fill, destroy, destroy_fill
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
                true, lat.size(), ident, fill, std::sqrt(2.)*matrix_elements[m], same_idx, pos1, create_count, create_fill_count, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                true, lat.size(), ident, fill, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, destroy_count, destroy_fill_count, create, create_fill
            ));

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, k, count, count));

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, j, e2d, d2e));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, matrix_elements[m], i, j, d2e, e2d));

            // here we have spin0--j--spin1--i--spin0
            // the sqrt(3.) counteracts the Clebsch coeff C^{110}_{mrm'} which applies when the spin1 couples back to spin0
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                false, lat.size(), ident, ident_full, std::sqrt(3.) * matrix_elements[m], i, j, flip_to_S0, flip_to_S2, flip_to_S0, flip_to_S2
            ));

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::two_term(false, ident, -0.5 * matrix_elements[m], i, j, count, count));

            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {

            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; k = i; l = j; }

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, std::sqrt(2.)*matrix_elements[m], same_idx, k, l, count, count_fill, create, create_fill, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, std::sqrt(2.)*matrix_elements[m], same_idx, l, k, count, count_fill, create, create_fill, destroy, destroy_fill
            ));

            used_elements[m] += 1;
        }

        // 9887 7371 8727

        // 4-fold degenerate (+spin) V_ijil = V_ijli = V_jiil = V_jili  <--- coded
        //                           V_ilij = V_ilji = V_liij = V_liji
        else if ( ((i==k && j!=l) || j==k || (j==l && i!=k)) && (i!=j && k!=l)) {
            int same_idx, pos1, pos2;
            if (i==k) { same_idx = i; pos1 = l; pos2 = j; }
            if (j==k) { same_idx = j; pos1 = l; pos2 = i; }
            if (j==l) { same_idx = j; pos1 = k; pos2 = i; }

            // Note: need minus because of clebsch gordan coeff from two destructors or two creators
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, e2d, e2d, destroy, destroy_fill, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, -std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, d2e, d2e, create, create_fill, create, create_fill
            ));

            if (! (same_idx < std::min(pos1,pos2))) continue;

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident_full, fill, std::sqrt(3./2.)*std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, flip_to_S2, flip_to_S2, create, create_fill_couple_down, destroy, destroy_fill
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos1, pos2, count, count, create, create_fill, destroy, destroy_fill
            ));

            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                // note minus sign, because commutation on same_idx is not taken into account
                ident_full, fill, -std::sqrt(3./2.)*std::sqrt(2.)*matrix_elements[m], same_idx, pos2, pos1, flip_to_S2, flip_to_S2, create, create_fill_couple_down, destroy, destroy_fill_couple_down
            ));
            this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::three_term(
                ident, fill, -0.5*std::sqrt(2.)*matrix_elements[m], same_idx, pos2, pos1, count, count, create, create_fill, destroy, destroy_fill
            ));

            used_elements[m] += 1;
        }

    } // matrix_elements for

    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
