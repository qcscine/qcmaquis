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

#include "dmrg/models/chem/term_maker.h"
#include "dmrg/models/chem/chem_detail.h"
#include "dmrg/models/chem/pg_util.h"

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
        return fill_cdagc;
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
    tag_type create_head, create_tail, destroy_head, destroy_tail,
             count, docc, e2d, d2e,
             ident, fill_cdagc, fill_ccdag, count_fill_cdagc, count_fill_ccdag;

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

    op_t ident_op;
    ident_op.insert_block(Matrix(1,1,1), A, A);
    ident_op.insert_block(Matrix(1,1,1), B, B);
    ident_op.insert_block(Matrix(1,1,1), C, C);
    ident_op.insert_block(Matrix(1,1,1), D, D);

    op_t fill_ccdag_op;
    fill_ccdag_op.insert_block(Matrix(1,1,1),  A, A);
    fill_ccdag_op.insert_block(Matrix(1,1,1),  D, D);
    fill_ccdag_op.insert_block(Matrix(1,1,-1), B, B);
    fill_ccdag_op.insert_block(Matrix(1,1,-1), C, C);
    fill_ccdag_op.insert_block(Matrix(1,1,1),  B, C);
    fill_ccdag_op.insert_block(Matrix(1,1,1),  C, B);

    op_t fill_cdagc_op;
    fill_cdagc_op.insert_block(Matrix(1,1,1),  A, A);
    fill_cdagc_op.insert_block(Matrix(1,1,1),  D, D);
    fill_cdagc_op.insert_block(Matrix(1,1,-1), B, B);
    fill_cdagc_op.insert_block(Matrix(1,1,-1), C, C);
    fill_cdagc_op.insert_block(Matrix(1,1,-1), B, C);
    fill_cdagc_op.insert_block(Matrix(1,1,-1), C, B);

    op_t create_tail_op;
    create_tail_op.insert_block(Matrix(1,1,-2.), B, A);
    create_tail_op.insert_block(Matrix(1,1,2.), C, A);
    create_tail_op.insert_block(Matrix(1,1,-sqrt(2.)), D, B);
    create_tail_op.insert_block(Matrix(1,1,sqrt(2.)), D, C);

    op_t destroy_tail_op;
    destroy_tail_op.insert_block(Matrix(1,1,1), A, B);
    destroy_tail_op.insert_block(Matrix(1,1,-1), A, C);
    destroy_tail_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
    destroy_tail_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

    op_t create_head_op;
    create_head_op.insert_block(Matrix(1,1,2.), B, A);
    create_head_op.insert_block(Matrix(1,1,2.), C, A);
    create_head_op.insert_block(Matrix(1,1,sqrt(2.)), D, B);
    create_head_op.insert_block(Matrix(1,1,sqrt(2.)), D, C);

    op_t destroy_head_op;
    destroy_head_op.insert_block(Matrix(1,1,1), A, B);
    destroy_head_op.insert_block(Matrix(1,1,1), A, C);
    destroy_head_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
    destroy_head_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

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

    op_t count_fill_cdagc_op;
    count_fill_cdagc_op.insert_block(Matrix(2,1,1),  A, A);
    count_fill_cdagc_op.insert_block(Matrix(1,1,-1), B, B);
    count_fill_cdagc_op.insert_block(Matrix(1,1,-1), C, C);
    count_fill_cdagc_op.insert_block(Matrix(1,1,-1), B, C);
    count_fill_cdagc_op.insert_block(Matrix(1,1,-1), C, B);

    op_t count_fill_ccdag_op;
    count_fill_ccdag_op.insert_block(Matrix(2,1,1),  A, A);
    count_fill_ccdag_op.insert_block(Matrix(1,1,-1), B, B);
    count_fill_ccdag_op.insert_block(Matrix(1,1,-1), C, C);
    count_fill_ccdag_op.insert_block(Matrix(1,1,1),  B, C);
    count_fill_ccdag_op.insert_block(Matrix(1,1,1),  C, B);

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);

    REGISTER(ident,        tag_detail::bosonic)
    REGISTER(fill_ccdag,   tag_detail::bosonic)
    REGISTER(fill_cdagc,   tag_detail::fermionic)
    REGISTER(create_head,  tag_detail::fermionic)
    REGISTER(create_tail,  tag_detail::fermionic)
    REGISTER(destroy_head, tag_detail::fermionic)
    REGISTER(destroy_tail, tag_detail::bosonic)
    REGISTER(count,        tag_detail::bosonic)
    REGISTER(docc,         tag_detail::bosonic)
    REGISTER(e2d,          tag_detail::bosonic)
    REGISTER(d2e,          tag_detail::bosonic)
    REGISTER(count_fill_cdagc,         tag_detail::bosonic)
    REGISTER(count_fill_ccdag,         tag_detail::bosonic)

#undef REGISTER
    /**********************************************************************/

#define PRINT(op) maquis::cout << #op << "\t" << op << std::endl;
    PRINT(ident)
    PRINT(fill_ccdag)
    PRINT(fill_cdagc)
    PRINT(create_head)
    PRINT(create_tail)
    PRINT(destroy_head)
    PRINT(destroy_tail)
    PRINT(count)
    PRINT(docc)
    PRINT(e2d)
    PRINT(d2e)
#undef PRINT

    chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, ident, ident, tag_handler);
    std::vector<value_type> & matrix_elements = term_assistant.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);
 
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = term_assistant.idx(m, 0);
        int j = term_assistant.idx(m, 1);
        int k = term_assistant.idx(m, 2);
        int l = term_assistant.idx(m, 3);

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
            int start = std::min(i,j), end = std::max(i,j);
            {
                term_descriptor term;
                term.is_fermionic = true;
                term.coeff = matrix_elements[m];

                for (int fs=0; fs < start; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(start, create_head) );

                for (int fs = start+1; fs < end; ++fs)
                    term.push_back( boost::make_tuple(fs, fill_cdagc) );
                term.push_back( boost::make_tuple(end, destroy_head) );

                for (int fs = end+1; fs < lat.size(); ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );

                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.is_fermionic = true;
                term.coeff = matrix_elements[m];

                for (int fs=0; fs < start; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(start, destroy_tail) );

                for (int fs = start+1; fs < end; ++fs)
                    term.push_back( boost::make_tuple(fs, fill_ccdag) );
                term.push_back( boost::make_tuple(end, create_tail) );

                for (int fs = end+1; fs < lat.size(); ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );

                this->terms_.push_back(term);
            }

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

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count, count);

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {

            term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, e2d, d2e);
            term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, d2e, e2d);
            //term_assistant.add_term(this->terms_,  -matrix_elements[m], i, j, count_single, count_single);

            // --> -c_j_up * cdag_j_down * c_i_down * cdag_i_up
            //term_assistant.add_term(
            //    this->terms_, matrix_elements[m], i, j, flip, flip
            //);

            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {

            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; k = i; l = j; }

            if ( same_idx > k || same_idx > l) continue;

            //int start = std::min(std::min(k,l), same_idx), end = std::max(std::max(k,l), same_idx);
            int start = same_idx, mid = std::min(k,l), end = std::max(k,l);
            {
                term_descriptor term;
                term.is_fermionic = true;
                term.coeff = matrix_elements[m];

                for (int fs=0; fs < start; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(start, count) );

                for (int fs = start+1; fs < mid; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(mid, create_head) );

                for (int fs = mid+1; fs < end; ++fs)
                    term.push_back( boost::make_tuple(fs, fill_cdagc) );
                term.push_back( boost::make_tuple(end, destroy_head) );

                for (int fs = end+1; fs < lat.size(); ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );

                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.is_fermionic = true;
                term.coeff = matrix_elements[m];

                for (int fs=0; fs < start; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(start, count) );

                for (int fs = start+1; fs < mid; ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );
                term.push_back( boost::make_tuple(mid, destroy_tail) );

                for (int fs = mid+1; fs < end; ++fs)
                    term.push_back( boost::make_tuple(fs, fill_ccdag) );
                term.push_back( boost::make_tuple(end, create_tail) );

                for (int fs = end+1; fs < lat.size(); ++fs)
                    term.push_back( boost::make_tuple(fs, ident) );

                this->terms_.push_back(term);
            }

            used_elements[m] += 1;
        }

    } // matrix_elements for

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
