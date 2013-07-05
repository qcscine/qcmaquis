/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2013 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_HAMILTONIANS_HPP
#define QC_HAMILTONIANS_HPP

template <class Matrix>
template <class M>
Hamiltonian<M, TwoU1> qc_model<Matrix>::H_impl() const
{
    typedef typename op_tag_t<M>::type tag_t;

    std::vector<typename hamterm_t<M>::type > terms;
    std::vector<int> used_elements(matrix_elements.size(), 0);

    typename op_t<M>::type create_up, create_down, destroy_up, destroy_down;
    typename op_t<M>::type count_up, count_down, docc, e2d, d2e;
    typename op_t<M>::type ident, fill;

    TwoU1::charge A(0), B(0), C(0), D(1);
    B[0]=1; C[1]=1;
    ident.insert_block(M(1, 1, 1), A, A);
    ident.insert_block(M(1, 1, 1), B, B);
    ident.insert_block(M(1, 1, 1), C, C);
    ident.insert_block(M(1, 1, 1), D, D);

    create_up.insert_block(M(1, 1, 1), A, B);
    create_up.insert_block(M(1, 1, 1), C, D);
    create_down.insert_block(M(1, 1, 1), A, C);
    create_down.insert_block(M(1, 1, 1), B, D);

    destroy_up.insert_block(M(1, 1, 1), B, A);
    destroy_up.insert_block(M(1, 1, 1), D, C);
    destroy_down.insert_block(M(1, 1, 1), C, A);
    destroy_down.insert_block(M(1, 1, 1), D, B);

    count_up.insert_block(M(1, 1, 1), B, B);
    count_up.insert_block(M(1, 1, 1), D, D);

    count_down.insert_block(M(1, 1, 1), C, C);
    count_down.insert_block(M(1, 1, 1), D, D);

    docc.insert_block(M(1, 1, 1), D, D);

    e2d.insert_block(M(1, 1, 1), A, D);
    d2e.insert_block(M(1, 1, 1), D, A);

    fill.insert_block(M(1, 1, 1), A, A);
    fill.insert_block(M(1, 1, -1), B, B);
    fill.insert_block(M(1, 1, -1), C, C);
    fill.insert_block(M(1, 1, 1), D, D);

    typename op_t<M>::type tmp;

    gemm(fill, create_down, tmp);
    create_down = tmp;
    gemm(destroy_down, fill, tmp);
    destroy_down = tmp;

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    boost::shared_ptr<OPTagTable<M, TwoU1> > op_table(new OPTagTable<M, TwoU1>());
    tag_t create_up_tag, create_down_tag, destroy_up_tag, destroy_down_tag,
          count_up_tag, count_down_tag, docc_tag, e2d_tag, d2e_tag,
          fill_tag, ident_tag;

    #define REGISTER(op) op ## _tag = op_table->register_op(op);

    REGISTER(ident)
    REGISTER(fill)
    REGISTER(create_up)
    REGISTER(create_down)
    REGISTER(destroy_up)
    REGISTER(destroy_down)
    REGISTER(count_up)
    REGISTER(count_down)
    REGISTER(e2d)
    REGISTER(d2e)
    REGISTER(docc)

    #undef REGISTER
    /**********************************************************************/
 
    /**********************************************************************/
    /*** Tagged terms *****************************************************/
    /**********************************************************************/
    std::vector<typename hamtagterm_t<M>::type > tagterms;
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = idx[m][0]-1;
        int j = idx[m][1]-1;
        int k = idx[m][2]-1;
        int l = idx[m][3]-1;

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {
            typename hamtagterm_t<M>::type term;
            term.fill_operator = ident_tag;
            term.scale = matrix_elements[m];
            term.operators.push_back(std::make_pair(0, ident_tag));
            tagterms.push_back(term);
            
            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                typename hamtagterm_t<M>::type term;
                term.fill_operator = ident_tag;
                term.scale = matrix_elements[m];
                term.operators.push_back(std::make_pair(i, count_up_tag));
                tagterms.push_back(term);
            }
            {
                typename hamtagterm_t<M>::type term;
                term.fill_operator = ident_tag;
                term.scale = matrix_elements[m];
                term.operators.push_back(std::make_pair(i, count_down_tag));
                tagterms.push_back(term);
            }

            used_elements[m] += 1;
            continue;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            tagterms.push_back(TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m], i, j, create_up_tag, destroy_up_tag, op_table));
            tagterms.push_back(TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m], i, j, create_down_tag, destroy_down_tag, op_table));
            tagterms.push_back(TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m], j, i, create_up_tag, destroy_up_tag, op_table));
            tagterms.push_back(TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m], j, i, create_down_tag, destroy_down_tag, op_table));

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {
            typename hamtagterm_t<M>::type term;
            term.fill_operator = ident_tag;
            term.scale = matrix_elements[m];
            term.operators.push_back(std::make_pair(i, docc_tag));
            tagterms.push_back(term);

            used_elements[m] += 1;
        }

        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            int same_idx, pos1;
            typename op_t<M>::type tmp;

            if ( i==j) { same_idx = i; pos1 = l; }
            if ( k==l) { same_idx = l; pos1 = i; }

            std::pair<tag_t, value_type> ptag;

            // 1a
            // --> c_l_up * n_i_down * cdag_i_up
            ptag = op_table->get_product_tag(count_down_tag, create_up_tag);
            tagterms.push_back( TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, destroy_up_tag, op_table) );

            // 1a_dagger
            // --> c_i_up * n_i_down * cdag_l_up
            ptag = op_table->get_product_tag(destroy_up_tag, count_down_tag);
            tagterms.push_back( TermMaker<M>::positional_two_term(true, fill_tag, -matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, create_up_tag, op_table) );

            // 1b
            // --> c_l_down * n_i_up * cdag_i_down (1b)
            ptag = op_table->get_product_tag(count_up_tag, create_down_tag);
            tagterms.push_back( TermMaker<M>::positional_two_term(true, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, destroy_down_tag, op_table) );

            // (1b)_dagger
            // --> c_i_down * n_i_up * cdag_l_down
            ptag = op_table->get_product_tag(destroy_down_tag, count_up_tag);
            tagterms.push_back( TermMaker<M>::positional_two_term(true, fill_tag, -matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, create_down_tag, op_table) );

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                maquis::cout << "matrix element underflow:\t" << std::abs(matrix_elements[m]) << std::endl;
                continue;
            }

            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, k,
                                            count_up_tag, count_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, k,
                                            count_up_tag, count_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, k,
                                            count_down_tag, count_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, k,
                                            count_down_tag, count_down_tag, op_table) );

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                maquis::cout << "matrix element underflow:\t" << std::abs(matrix_elements[m]) << std::endl;
                used_elements[m] += 1;
                continue;
            }
            typename op_t<M>::type tmp1, tmp2;

            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, j,
                                            e2d_tag, d2e_tag, op_table) );
            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, matrix_elements[m], i, j,
                                            d2e_tag, e2d_tag, op_table) );

            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, -matrix_elements[m], i, j,
                                            count_up_tag, count_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::two_term(false, ident_tag, -matrix_elements[m], i, j,
                                            count_down_tag, count_down_tag, op_table) );

            std::pair<tag_t, value_type> ptag1, ptag2;

            // Could insert fill operators without changing the result
            // --> -c_j_up * cdag_j_down * c_i_down * cdag_i_up
            ptag1 = op_table->get_product_tag(destroy_down_tag, create_up_tag);
            ptag2 = op_table->get_product_tag(destroy_up_tag, create_down_tag);
            tagterms.push_back(TermMaker<M>::two_term(false, ident_tag, -matrix_elements[m] * ptag1.second * ptag2.second, i, j,
                            ptag1.first, ptag2.first, op_table));

            // --> -c_i_up * cdag_i_down * c_j_down * cdag_j_up
            ptag1 = op_table->get_product_tag(destroy_up_tag, create_down_tag);
            ptag2 = op_table->get_product_tag(destroy_down_tag, create_up_tag);
            tagterms.push_back(TermMaker<M>::two_term(false, ident_tag, -matrix_elements[m] * ptag1.second * ptag2.second, i, j,
                            ptag1.first, ptag2.first, op_table));
            
            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; k = i; l = j; }

            // n_up * cdag_up * c_up <--
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, k, l,
                                             create_up_tag, destroy_up_tag, create_up_tag, destroy_up_tag, op_table));
            // n_up * cdag_down * c_down <--
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, k, l,
                                             create_up_tag, destroy_up_tag, create_down_tag, destroy_down_tag, op_table));
            // n_down * cdag_up * c_up <--
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, k, l,
                                             create_down_tag, destroy_down_tag, create_up_tag, destroy_up_tag, op_table));
            // n_down * cdag_down * c_down <--
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, k, l,
                                             create_down_tag, destroy_down_tag, create_down_tag, destroy_down_tag, op_table));

            // --> n_up * c_up * cdag_up
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, l, k,
                                             create_up_tag, destroy_up_tag, create_up_tag, destroy_up_tag, op_table));
            // --> n_up * c_down * cdag_down
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, l, k,
                                             create_up_tag, destroy_up_tag, create_down_tag, destroy_down_tag, op_table));
            // --> n_down * c_up * cdag_up
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, l, k,
                                             create_down_tag, destroy_down_tag, create_up_tag, destroy_up_tag, op_table));
            // --> n_down * c_down * cdag_down
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m], same_idx, l, k,
                                             create_down_tag, destroy_down_tag, create_down_tag, destroy_down_tag, op_table));

            used_elements[m] += 1;
        }

        // 9887 7371 8727

        // 4-fold degenerate (+spin) V_ijil = V_ijli = V_jiil = V_jili  <--- coded
        //                           V_ilij = V_ilji = V_liij = V_liji
        else if ( ((i==k && j!=l) || j==k || (j==l && i!=k)) && (i!=j && k!=l)) {
        //else if ( i==2 && (((i==k && j!=l) || j==k || (j==l && i!=k)) && (i!=j && k!=l)) ) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            int same_idx, pos1, pos2;
            if (i==k) { same_idx = i; pos1 = l; pos2 = j; }
            if (j==k) { same_idx = j; pos1 = l; pos2 = i; }
            if (j==l) { same_idx = j; pos1 = k; pos2 = i; }

            std::pair<tag_t, value_type> ptag;

            ptag = op_table->get_product_tag(create_up_tag, fill_tag);
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1, pos2,
                                             ptag.first, create_down_tag, destroy_down_tag, destroy_up_tag, op_table) );
            ptag = op_table->get_product_tag(create_down_tag, fill_tag);
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1, pos2,
                                             ptag.first, create_up_tag, destroy_up_tag, destroy_down_tag, op_table) );
            ptag = op_table->get_product_tag(destroy_down_tag, fill_tag);
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1, pos2,
                                             ptag.first, destroy_up_tag, create_up_tag, create_down_tag, op_table) );
            ptag = op_table->get_product_tag(destroy_up_tag, fill_tag);
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, matrix_elements[m] * ptag.second, same_idx, pos1, pos2,
                                             ptag.first, destroy_down_tag, create_down_tag, create_up_tag, op_table) );

            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos1, pos2,
                                             create_up_tag, destroy_up_tag, create_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos1, pos2,
                                             create_up_tag, destroy_down_tag, create_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos1, pos2,
                                             create_down_tag, destroy_up_tag, create_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos1, pos2,
                                             create_down_tag, destroy_down_tag, create_down_tag, destroy_down_tag, op_table) );

            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos2, pos1,
                                             create_up_tag, destroy_up_tag, create_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos2, pos1,
                                             create_up_tag, destroy_down_tag, create_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos2, pos1,
                                             create_down_tag, destroy_up_tag, create_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::three_term(true, ident_tag, fill_tag, -matrix_elements[m], same_idx, pos2, pos1,
                                             create_down_tag, destroy_down_tag, create_down_tag, destroy_down_tag, op_table) );

            used_elements[m] += 1;
        }

        // 32 (8x4)-fold degenerate V_ijkl = V_jikl = V_ijlk = V_jilk = V_klij = V_lkij = V_klji = V_lkji * spin
        // V_ijkl -> 24 permutations which fall into 3 equivalence classes of 8 elements (with identical V_ matrix element)
        // coded: 4 index permutations x 4 spin combinations 
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {
        //else if ( i==4 && i!=j && j!=k && k!=l && i!=k && j!=l) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            // 1
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,k,l,j,
                                create_up_tag, create_up_tag, destroy_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,k,l,j,
                                create_up_tag, create_down_tag, destroy_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,k,l,j,
                                create_down_tag, create_up_tag, destroy_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,k,l,j,
                                create_down_tag, create_down_tag, destroy_down_tag, destroy_down_tag, op_table) );

            // 2
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,l,k,j,
                                create_up_tag, create_up_tag, destroy_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,l,k,j,
                                create_up_tag, create_down_tag, destroy_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,l,k,j,
                                create_down_tag, create_up_tag, destroy_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], i,l,k,j,
                                create_down_tag, create_down_tag, destroy_down_tag, destroy_down_tag, op_table) );

            // 3
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,k,l,i,
                                create_up_tag, create_up_tag, destroy_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,k,l,i,
                                create_up_tag, create_down_tag, destroy_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,k,l,i,
                                create_down_tag, create_up_tag, destroy_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,k,l,i,
                                create_down_tag, create_down_tag, destroy_down_tag, destroy_down_tag, op_table) );

            // 4
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,l,k,i,
                                create_up_tag, create_up_tag, destroy_up_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,l,k,i,
                                create_up_tag, create_down_tag, destroy_down_tag, destroy_up_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,l,k,i,
                                create_down_tag, create_up_tag, destroy_up_tag, destroy_down_tag, op_table) );
            tagterms.push_back( TermMaker<M>::four_term(true, ident_tag, fill_tag, matrix_elements[m], j,l,k,i,
                                create_down_tag, create_down_tag, destroy_down_tag, destroy_down_tag, op_table) );

            used_elements[m] += 1;
        }
    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    assert( it_0 == used_elements.end() );

    maquis::cout << "The hamiltonian will contain " << tagterms.size() << " terms\n";

    return Hamiltonian<M, TwoU1>(phys, ident, terms, ident_tag, tagterms, op_table);

}
    
#endif
