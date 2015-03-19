/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2015 Institute for Theoretical Physics, ETH Zurich
 *               2012-2015 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_MODEL_HPP
#define QC_MODEL_HPP

template <class Matrix, class SymmGroup>
qc_model<Matrix, SymmGroup>::qc_model(Lattice const & lat_, BaseParameters & parms_)
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

    typename SymmGroup::charge A(0), B(0), C(0), D(1);
    B[0]=1; C[1]=1;

    for (subcharge irr=0; irr <= max_irrep; ++irr)
    {
        Index<SymmGroup> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(B, irr), 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(C, irr), 1));
        phys.insert(std::make_pair(D, 1));

        phys_indices.push_back(phys);
    }

    op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
         count_up_op, count_down_op, count_up_down_op, docc_op, e2d_op, d2e_op,
         ident_op, fill_op;

    ident_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_op.insert_block(Matrix(1, 1, 1), B, B);
    ident_op.insert_block(Matrix(1, 1, 1), C, C);
    ident_op.insert_block(Matrix(1, 1, 1), D, D);

    create_up_op.insert_block(Matrix(1, 1, 1), A, B);
    create_up_op.insert_block(Matrix(1, 1, 1), C, D);
    create_down_op.insert_block(Matrix(1, 1, 1), A, C);
    create_down_op.insert_block(Matrix(1, 1, 1), B, D);

    destroy_up_op.insert_block(Matrix(1, 1, 1), B, A);
    destroy_up_op.insert_block(Matrix(1, 1, 1), D, C);
    destroy_down_op.insert_block(Matrix(1, 1, 1), C, A);
    destroy_down_op.insert_block(Matrix(1, 1, 1), D, B);

    count_up_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_op.insert_block(Matrix(1, 1, 1), D, D);

    count_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_down_op.insert_block(Matrix(1, 1, 1), D, D);

    count_up_down_op.insert_block(Matrix(1, 1, 1), B, B);
    count_up_down_op.insert_block(Matrix(1, 1, 1), C, C);
    count_up_down_op.insert_block(Matrix(1, 1, 2), D, D);

    docc_op.insert_block(Matrix(1, 1, 1), D, D);

    e2d_op.insert_block(Matrix(1, 1, 1), A, D);
    d2e_op.insert_block(Matrix(1, 1, 1), D, A);

    fill_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_op.insert_block(Matrix(1, 1, -1), B, B);
    fill_op.insert_block(Matrix(1, 1, -1), C, C);
    fill_op.insert_block(Matrix(1, 1, 1), D, D);

    op_t tmp;

    gemm(fill_op, create_down_op, tmp);
    create_down_op = tmp;
    gemm(destroy_down_op, fill_op, tmp);
    destroy_down_op = tmp;

    // only effective if point group symmetry is active, need to adapt operators to different irreps
    #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

    GENERATE_SITE_SPECIFIC(ident_op)
    GENERATE_SITE_SPECIFIC(fill_op)
    GENERATE_SITE_SPECIFIC(create_up_op)
    GENERATE_SITE_SPECIFIC(create_down_op)
    GENERATE_SITE_SPECIFIC(destroy_up_op)
    GENERATE_SITE_SPECIFIC(destroy_down_op)
    GENERATE_SITE_SPECIFIC(count_up_op)
    GENERATE_SITE_SPECIFIC(count_down_op)
    GENERATE_SITE_SPECIFIC(e2d_op)
    GENERATE_SITE_SPECIFIC(d2e_op)
    GENERATE_SITE_SPECIFIC(docc_op)

    #undef GENERATE_SITE_SPECIFIC

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    #define REGISTER(op, kind) op = this->register_site_specific(op ## _ops, kind);

    REGISTER(ident,        tag_detail::bosonic)
    REGISTER(fill,         tag_detail::bosonic)
    REGISTER(create_up,    tag_detail::fermionic)
    REGISTER(create_down,  tag_detail::fermionic)
    REGISTER(destroy_up,   tag_detail::fermionic)
    REGISTER(destroy_down, tag_detail::fermionic)
    REGISTER(count_up,     tag_detail::bosonic)
    REGISTER(count_down,   tag_detail::bosonic)
    REGISTER(e2d,          tag_detail::bosonic)
    REGISTER(d2e,          tag_detail::bosonic)
    REGISTER(docc,         tag_detail::bosonic)
    REGISTER(count_up_down,         tag_detail::bosonic)

    #undef REGISTER
    /**********************************************************************/
    std::pair<std::vector<tag_type>, std::vector<value_type> > create_up_times_fill = tag_handler->get_product_tags(create_up, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > create_down_times_fill = tag_handler->get_product_tags(create_down, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > destroy_up_times_fill = tag_handler->get_product_tags(destroy_up, fill);
    std::pair<std::vector<tag_type>, std::vector<value_type> > destroy_down_times_fill = tag_handler->get_product_tags(destroy_down, fill);
    /**********************************************************************/

    chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, ident, fill, tag_handler);
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
            term.push_back( boost::make_tuple(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]));
            this->terms_.push_back(term);
            
            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                term_descriptor term;
                term.coeff = matrix_elements[m];
                term.push_back( boost::make_tuple(i, count_up[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.coeff = matrix_elements[m];
                term.push_back( boost::make_tuple(i, count_down[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
                this->terms_.push_back(term);
            }

            used_elements[m] += 1;
            continue;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {

            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill, matrix_elements[m], i, j, create_up, destroy_up, tag_handler, lat)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill, matrix_elements[m], i, j, create_down, destroy_down, tag_handler, lat)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill, matrix_elements[m], j, i, create_up, destroy_up, tag_handler, lat)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill, matrix_elements[m], j, i, create_down, destroy_down, tag_handler, lat)
            );

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back(boost::make_tuple(i, docc[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]));
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {

            int same_idx, pos1;

            if      (i==j) { same_idx = i; pos1 = l; }
            else if (k==l) { same_idx = l; pos1 = i; }
            else           { throw std::runtime_error("Term generation logic has failed for V_ijjj term\n"); }

            //std::pair<std::vector<tag_type>, value_type> ptag;

            // 1a
            // --> c_l_up * n_i_down * cdag_i_up
            //ptag = tag_handler->get_product_tags(count_down, create_up);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, matrix_elements[m], same_idx, pos1,
                                           count_down, create_up, destroy_up, tag_handler, lat) );

            // 1a_dagger
            // --> c_i_up * n_i_down * cdag_l_up
            //ptag = tag_handler->get_product_tags(destroy_up, count_down);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, -matrix_elements[m], same_idx, pos1,
                                           destroy_up, count_down, create_up, tag_handler, lat) );

            // 1b
            // --> c_l_down * n_i_up * cdag_i_down (1b)
            //ptag = tag_handler->get_product_tags(count_up, create_down);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, matrix_elements[m], same_idx, pos1,
                                           count_up, create_down, destroy_down, tag_handler, lat) );

            // (1b)_dagger
            // --> c_i_down * n_i_up * cdag_l_down
            //ptag = tag_handler->get_product_tags(destroy_down, count_up);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, -matrix_elements[m], same_idx, pos1,
                                         destroy_down, count_up , create_down, tag_handler, lat) );

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            //term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_up, count_up);
            //term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_up, count_down);
            //term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_down, count_up);
            //term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_down, count_down);
            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_up_down, count_up_down);

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {

            term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, e2d, d2e);
            term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, d2e, e2d);
            term_assistant.add_term(this->terms_, -matrix_elements[m], i, j, count_up, count_up);
            term_assistant.add_term(this->terms_, -matrix_elements[m], i, j, count_down, count_down);

            std::pair<std::vector<tag_type>, value_type> ptag1, ptag2;

            // Could insert fill operators without changing the result
            // --> -c_j_up * cdag_j_down * c_i_down * cdag_i_up
            //ptag1 = tag_handler->get_product_tags(destroy_down, create_up);
            //ptag2 = tag_handler->get_product_tags(destroy_up, create_down);
            term_assistant.add_term(
                //this->terms_, -matrix_elements[m] * ptag1.second * ptag2.second, i, j, ptag1.first, ptag2.first
                this->terms_, -matrix_elements[m], i, j, destroy_down, create_up, destroy_up, create_down
            );

            // --> -c_i_up * cdag_i_down * c_j_down * cdag_j_up
            //ptag1 = tag_handler->get_product_tags(destroy_up, create_down);
            //ptag2 = tag_handler->get_product_tags(destroy_down, create_up);
            term_assistant.add_term(
                //this->terms_, -matrix_elements[m] * ptag1.second * ptag2.second, i, j, ptag1.first, ptag2.first
                this->terms_, -matrix_elements[m], i, j, destroy_up, create_down, destroy_down, create_up
            );
            
            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {

            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; k = i; l = j; }

            // n_up * cdag_up * c_up <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_up, destroy_up, create_up, destroy_up);
            // n_up * cdag_down * c_down <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_up, destroy_up, create_down, destroy_down);
            // n_down * cdag_up * c_up <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_down, destroy_down, create_up, destroy_up);
            // n_down * cdag_down * c_down <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_down, destroy_down, create_down, destroy_down);

            // --> n_up * c_up * cdag_up
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_up, destroy_up, create_up, destroy_up);
            // --> n_up * c_down * cdag_down
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_up, destroy_up, create_down, destroy_down);
            // --> n_down * c_up * cdag_up
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_down, destroy_down, create_up, destroy_up);
            // --> n_down * c_down * cdag_down
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_down, destroy_down, create_down, destroy_down);

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

            //std::pair<std::vector<tag_type>, value_type> ptag;
            typename SymmGroup::subcharge irr = lat.get_prop<typename SymmGroup::subcharge>("type", same_idx);

            //ptag = tag_handler->get_product_tags(create_up, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*create_up_times_fill.second[irr], same_idx, pos1, pos2, create_up_times_fill.first, create_down , destroy_down, destroy_up
                //this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, create_down , destroy_down, destroy_up
            );
            //ptag = tag_handler->get_product_tags(create_down, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*create_down_times_fill.second[irr], same_idx, pos1, pos2, create_down_times_fill.first, create_up   , destroy_up  , destroy_down
                //this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, create_up   , destroy_up  , destroy_down
            );
            //ptag = tag_handler->get_product_tags(destroy_down, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*destroy_down_times_fill.second[irr], same_idx, pos1, pos2, destroy_down_times_fill.first, destroy_up  , create_up   , create_down
                //this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, destroy_up  , create_up   , create_down
            );
            //ptag = tag_handler->get_product_tags(destroy_up, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*destroy_up_times_fill.second[irr], same_idx, pos1, pos2, destroy_up_times_fill.first, destroy_down, create_down , create_up
                //this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, destroy_down, create_down , create_up
            );

            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_up,   destroy_up,   create_up,   destroy_up
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_up,   destroy_down, create_down, destroy_up
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_down, destroy_up,   create_up,   destroy_down
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_down, destroy_down, create_down, destroy_down
            );

            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_up,   destroy_up,   create_up,   destroy_up
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_up,   destroy_down, create_down, destroy_up
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_down, destroy_up,   create_up,   destroy_down
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_down, destroy_down, create_down, destroy_down
            );

            used_elements[m] += 1;
        }

        // 32 (8x4)-fold degenerate V_ijkl = V_jikl = V_ijlk = V_jilk = V_klij = V_lkij = V_klji = V_lkji * spin
        // V_ijkl -> 24 permutations which fall into 3 equivalence classes of 8 elements (with identical V_ matrix element)
        // coded: 4 index permutations x 4 spin combinations 
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {
            
            // 1
            term_assistant.add_term(this->terms_, i,k,l,j, create_up, create_up, destroy_up, destroy_up);
            term_assistant.add_term(this->terms_, i,k,l,j, create_up, create_down, destroy_down, destroy_up);
            term_assistant.add_term(this->terms_, i,k,l,j, create_down, create_up, destroy_up, destroy_down);
            term_assistant.add_term(this->terms_, i,k,l,j, create_down, create_down, destroy_down, destroy_down);

            // 2
            term_assistant.add_term(this->terms_, i,l,k,j, create_up, create_up, destroy_up, destroy_up);
            term_assistant.add_term(this->terms_, i,l,k,j, create_up, create_down, destroy_down, destroy_up);
            term_assistant.add_term(this->terms_, i,l,k,j, create_down, create_up, destroy_up, destroy_down);
            term_assistant.add_term(this->terms_, i,l,k,j, create_down, create_down, destroy_down, destroy_down);

            // 3
            term_assistant.add_term(this->terms_, j,k,l,i, create_up, create_up, destroy_up, destroy_up);
            term_assistant.add_term(this->terms_, j,k,l,i, create_up, create_down, destroy_down, destroy_up);
            term_assistant.add_term(this->terms_, j,k,l,i, create_down, create_up, destroy_up, destroy_down);
            term_assistant.add_term(this->terms_, j,k,l,i, create_down, create_down, destroy_down, destroy_down);

            // 4
            term_assistant.add_term(this->terms_, j,l,k,i, create_up, create_up, destroy_up, destroy_up);
            term_assistant.add_term(this->terms_, j,l,k,i, create_up, create_down, destroy_down, destroy_up);
            term_assistant.add_term(this->terms_, j,l,k,i, create_down, create_up, destroy_up, destroy_down);
            term_assistant.add_term(this->terms_, j,l,k,i, create_down, create_down, destroy_down, destroy_down);
        
            used_elements[m] += 1;
        }
    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    assert( it_0 == used_elements.end() );

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}
    
#endif
