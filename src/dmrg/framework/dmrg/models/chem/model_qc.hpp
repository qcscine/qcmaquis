/*****************************************************************************
 *
 * MAQUIS DMRG Project
 *
 * Copyright (C) 2012-2012 by Sebastian Keller <sebkelle@phys.ethz.ch>
 *
 *
 *****************************************************************************/

#ifndef QC_HAMILTONIANS_HPP
#define QC_HAMILTONIANS_HPP

template<class Matrix>
Hamiltonian<Matrix, TwoU1> qc_model<Matrix>::H() const
{
    std::vector<hamterm_t> terms;
    std::vector<int> used_elements(matrix_elements.size(), 0);
    
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = idx[m][0]-1;
        int j = idx[m][1]-1;
        int k = idx[m][2]-1;
        int l = idx[m][3]-1;

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {
            hamterm_t term;
            term.fill_operator = ident;
            term.operators.push_back( std::make_pair(0, matrix_elements[m]*ident) );
            terms.push_back(term);
            
            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                hamterm_t term;
                term.fill_operator = ident;
                term.operators.push_back( std::make_pair(i, matrix_elements[m]*count_up) );
                terms.push_back(term);
            }
            {
                hamterm_t term;
                term.fill_operator = ident;
                term.operators.push_back( std::make_pair(i, matrix_elements[m]*count_down) );
                terms.push_back(term);
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
            terms.push_back( make_positional_two_term(true, fill, i, j, matrix_elements[m]*create_up, destroy_up) );
            terms.push_back( make_positional_two_term(true, fill, i, j, matrix_elements[m]*create_down, destroy_down) );

            terms.push_back( make_positional_two_term(true, fill, j, i, matrix_elements[m]*create_up, destroy_up) );
            terms.push_back( make_positional_two_term(true, fill, j, i, matrix_elements[m]*create_down, destroy_down) );

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {
            hamterm_t term;
            term.fill_operator = ident;
            term.operators.push_back( std::make_pair(i, matrix_elements[m]*doubly_occ) );
            terms.push_back(term);

            used_elements[m] += 1;
        }

        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                continue;
            }
            int same_idx, pos1;
            op_t tmp;

            if ( i==j) { same_idx = i; pos1 = l; }
            if ( k==l) { same_idx = l; pos1 = i; }

            // 1a
            // --> c_l_up * n_i_down * cdag_i_up
            gemm(count_down, create_up, tmp);
            terms.push_back( make_positional_two_term(true, fill, same_idx, pos1,
                                           matrix_elements[m]*tmp, destroy_up) );
            // 1a_dagger
            // --> c_i_up * n_i_down * cdag_l_up
            gemm(destroy_up, count_down, tmp);
            terms.push_back( make_positional_two_term(true, fill, same_idx, pos1,
                                           -1*matrix_elements[m]*tmp, create_up) );

            // 1b
            // --> c_l_down * n_i_up * cdag_i_down (1b)
            gemm(count_up, create_down, tmp);
            terms.push_back( make_positional_two_term(true, fill, same_idx, pos1,
                                           matrix_elements[m]*tmp, destroy_down) );

            // (1b)_dagger
            // --> c_i_down * n_i_up * cdag_l_down
            gemm(destroy_down, count_up, tmp);
            terms.push_back( make_positional_two_term(true, fill, same_idx, pos1,
                                           -1*matrix_elements[m]*tmp, create_down) );

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                used_elements[m] += 1;
                maquis::cout << "matrix element underflow:\t" << std::abs(matrix_elements[m]) << std::endl;
                continue;
            }

            terms.push_back( make_two_term(false, ident, i, k,
                                            matrix_elements[m]*count_up, count_up) );
            terms.push_back( make_two_term(false, ident, i, k,
                                            matrix_elements[m]*count_up, count_down) );
            terms.push_back( make_two_term(false, ident, i, k,
                                            matrix_elements[m]*count_down, count_up) );
            terms.push_back( make_two_term(false, ident, i, k,
                                            matrix_elements[m]*count_down, count_down) );

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {
            if (std::abs(matrix_elements[m]) < 1.e-20) {
                maquis::cout << "matrix element underflow:\t" << std::abs(matrix_elements[m]) << std::endl;
                used_elements[m] += 1;
                continue;
            }
            op_t tmp1, tmp2;

            terms.push_back( make_two_term(false, ident, i, j,
                                            matrix_elements[m]*empty2doubly_occ, doubly_occ2empty) );
            terms.push_back( make_two_term(false, ident, i, j,
                                            matrix_elements[m]*doubly_occ2empty, empty2doubly_occ) );

            terms.push_back( make_two_term(false, ident, i, j,
                                            -1*matrix_elements[m]*count_up, count_up) );
            terms.push_back( make_two_term(false, ident, i, j,
                                            -1*matrix_elements[m]*count_down, count_down) );

            // Could insert fill operators without changing the result
            // --> -c_j_up * cdag_j_down * c_i_down * cdag_i_up
            gemm(destroy_down, create_up, tmp1);
            gemm(destroy_up, create_down, tmp2);
            terms.push_back( make_two_term(false, ident, i, j,
                                            -1*matrix_elements[m]*tmp1, tmp2) );

            // --> -c_i_up * cdag_i_down * c_j_down * cdag_j_up
            gemm(destroy_up, create_down, tmp1);
            gemm(destroy_down, create_up, tmp2);
            terms.push_back( make_two_term(false, ident, i, j,
                                            -1*matrix_elements[m]*tmp1, tmp2) );
            
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
            terms.push_back( make_three_term(true, fill, same_idx, k, l,
                                             matrix_elements[m]*create_up, destroy_up, create_up, destroy_up));
            // n_up * cdag_down * c_down <--
            terms.push_back( make_three_term(true, fill, same_idx, k, l,
                                             matrix_elements[m]*create_up, destroy_up, create_down, destroy_down));
            // n_down * cdag_up * c_up <--
            terms.push_back( make_three_term(true, fill, same_idx, k, l,
                                             matrix_elements[m]*create_down, destroy_down, create_up, destroy_up));
            // n_down * cdag_down * c_down <--
            terms.push_back( make_three_term(true, fill, same_idx, k, l,
                                             matrix_elements[m]*create_down, destroy_down, create_down, destroy_down));

            // --> n_up * c_up * cdag_up
            terms.push_back( make_three_term(true, fill, same_idx, l, k,
                                             matrix_elements[m]*create_up, destroy_up, create_up, destroy_up));
            // --> n_up * c_down * cdag_down
            terms.push_back( make_three_term(true, fill, same_idx, l, k,
                                             matrix_elements[m]*create_up, destroy_up, create_down, destroy_down));
            // --> n_down * c_up * cdag_up
            terms.push_back( make_three_term(true, fill, same_idx, l, k,
                                             matrix_elements[m]*create_down, destroy_down, create_up, destroy_up));
            // --> n_down * c_down * cdag_down
            terms.push_back( make_three_term(true, fill, same_idx, l, k,
                                             matrix_elements[m]*create_down, destroy_down, create_down, destroy_down));

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
            op_t tmp;

            gemm(create_up, fill, tmp);
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             matrix_elements[m]*tmp, create_down, destroy_down, destroy_up) );
            gemm(create_down, fill, tmp);
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             matrix_elements[m]*tmp, create_up, destroy_up, destroy_down) );
            gemm(destroy_down, fill, tmp);
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             matrix_elements[m]*tmp, destroy_up, create_up, create_down) );
            gemm(destroy_up, fill, tmp);
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             matrix_elements[m]*tmp, destroy_down, create_down, create_up) );

            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             -1*matrix_elements[m]*create_up, destroy_up, create_up, destroy_up) );
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             -1*matrix_elements[m]*create_up, destroy_down, create_down, destroy_up) );
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             -1*matrix_elements[m]*create_down, destroy_up, create_up, destroy_down) );
            terms.push_back( make_three_term(true, fill, same_idx, pos1, pos2,
                                             -1*matrix_elements[m]*create_down, destroy_down, create_down, destroy_down) );

            terms.push_back( make_three_term(true, fill, same_idx, pos2, pos1,
                                             -1*matrix_elements[m]*create_up, destroy_up, create_up, destroy_up) );
            terms.push_back( make_three_term(true, fill, same_idx, pos2, pos1,
                                             -1*matrix_elements[m]*create_up, destroy_down, create_down, destroy_up) );
            terms.push_back( make_three_term(true, fill, same_idx, pos2, pos1,
                                             -1*matrix_elements[m]*create_down, destroy_up, create_up, destroy_down) );
            terms.push_back( make_three_term(true, fill, same_idx, pos2, pos1,
                                             -1*matrix_elements[m]*create_down, destroy_down, create_down, destroy_down) );

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
            terms.push_back( make_positional_four_term(true, fill, i,k,l,j,
                             matrix_elements[m]*create_up, create_up, destroy_up, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, i,k,l,j,
                             matrix_elements[m]*create_up, create_down, destroy_down, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, i,k,l,j,
                             matrix_elements[m]*create_down, create_up, destroy_up, destroy_down) );
            terms.push_back( make_positional_four_term(true, fill, i,k,l,j,
                             matrix_elements[m]*create_down, create_down, destroy_down, destroy_down) );

            // 2
            terms.push_back( make_positional_four_term(true, fill, i,l,k,j,
                             matrix_elements[m]*create_up, create_up, destroy_up, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, i,l,k,j,
                             matrix_elements[m]*create_up, create_down, destroy_down, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, i,l,k,j,
                             matrix_elements[m]*create_down, create_up, destroy_up, destroy_down) );
            terms.push_back( make_positional_four_term(true, fill, i,l,k,j,
                             matrix_elements[m]*create_down, create_down, destroy_down, destroy_down) );

            // 3
            terms.push_back( make_positional_four_term(true, fill, j,k,l,i,
                             matrix_elements[m]*create_up, create_up, destroy_up, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, j,k,l,i,
                             matrix_elements[m]*create_up, create_down, destroy_down, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, j,k,l,i,
                             matrix_elements[m]*create_down, create_up, destroy_up, destroy_down) );
            terms.push_back( make_positional_four_term(true, fill, j,k,l,i,
                             matrix_elements[m]*create_down, create_down, destroy_down, destroy_down) );

            // 4
            terms.push_back( make_positional_four_term(true, fill, j,l,k,i,
                             matrix_elements[m]*create_up, create_up, destroy_up, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, j,l,k,i,
                             matrix_elements[m]*create_up, create_down, destroy_down, destroy_up) );
            terms.push_back( make_positional_four_term(true, fill, j,l,k,i,
                             matrix_elements[m]*create_down, create_up, destroy_up, destroy_down) );
            terms.push_back( make_positional_four_term(true, fill, j,l,k,i,
                             matrix_elements[m]*create_down, create_down, destroy_down, destroy_down) );

            used_elements[m] += 1;
        }
    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    assert( it_0 == used_elements.end() );

    maquis::cout << "The hamiltonian will contain " << terms.size() << " terms\n";

    return ham(phys, ident, terms);

}
    
#endif
