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

#ifndef REL_QC_MODEL_HPP
#define REL_QC_MODEL_HPP

template <class Matrix, class SymmGroup>
//qc_model<Matrix, SymmGroup>::qc_model(Lattice const & lat_, BaseParameters & parms_)
rel_qc_model<Matrix, SymmGroup>::rel_qc_model(Lattice const & lat_, BaseParameters & parms_)
: lat(lat_)
, parms(parms_)
, tag_handler(new table_type())
{
    // find the highest irreducible representation number 
    // used to generate ops for all irreps 0..max_irrep
    max_irrep = 0;
    for (pos_t p=0; p < lat.size(); ++p)
        max_irrep = (lat.get_prop<typename SymmGroup::subcharge>("type", p) > max_irrep)
                        ? lat.get_prop<typename SymmGroup::subcharge>("type", p) : max_irrep;

    typename SymmGroup::charge A(0), B(0), C(0), D(1);
    B[0]=1; C[1]=1;

    // A -> empty, B -> unbarred spinor, C-> barred spinor
    phys.insert(std::make_pair(A, 1));
    phys.insert(std::make_pair(B, 1));
    phys.insert(std::make_pair(C, 1));

    // create physical indices
    phys_indices.resize(lat.size());
    for (std::size_t site = 0; site < lat.size(); ++site)
    {
        // Set point group
        int pgrp = lat.get_prop<int>("irrep", site);

        B[2] = pgrp;
        C[2] = pgrp;

        Index<SymmGroup> loc;
        if (site < lat.size()/2)
        {   // create unbarred bases in the first half
            loc.insert(std::make_pair(A, 1));
            loc.insert(std::make_pair(B, 1));
        }
        else
        {   // barred bases in the second half
            loc.insert(std::make_pair(A, 1));
            loc.insert(std::make_pair(C, 1));
        }
        phys_indices[site] = loc;
    }

    op_t create_unbar_op, create_bar_op, destroy_unbar_op, destroy_bar_op,
         count_unbar_op, count_bar_op, docc_op, e2d_op, d2e_op,
         ident_unbar_op, ident_bar_op, fill_unbar_op, fill_bar_op;

    ident_unbar_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_unbar_op.insert_block(Matrix(1, 1, 1), B, B);
    
    ident_bar_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_bar_op.insert_block(Matrix(1, 1, 1), C, C);

    create_unbar_op.insert_block(Matrix(1, 1, 1), A, B);
    create_bar_op.insert_block(Matrix(1, 1, 1), A, C);

    destroy_unbar_op.insert_block(Matrix(1, 1, 1), B, A);
    destroy_bar_op.insert_block(Matrix(1, 1, 1), C, A);

    count_unbar_op.insert_block(Matrix(1, 1, 1), B, B);
    count_bar_op.insert_block(Matrix(1, 1, 1), C, C);
 
    docc_op.insert_block(Matrix(1, 1, 1), B, B);
    docc_op.insert_block(Matrix(1, 1, 1), C, C);
    
    fill_unbar_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_unbar_op.insert_block(Matrix(1, 1, -1), B, B);
    
    fill_bar_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_bar_op.insert_block(Matrix(1, 1, -1), C, C);

    op_t tmp;
    tag_type dummy;

    // Not sure if fill_bar is right
    gemm(fill_bar_op, create_bar_op, tmp);
    create_bar_op = tmp;
    gemm(destroy_bar_op, fill_bar_op, tmp);
    destroy_bar_op = tmp;

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    #define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);

    REGISTER(ident_unbar,    tag_detail::bosonic)
    REGISTER(ident_bar,      tag_detail::bosonic)
    REGISTER(fill_unbar,     tag_detail::bosonic)
    REGISTER(fill_bar,       tag_detail::bosonic)
    REGISTER(create_unbar,   tag_detail::fermionic)
    REGISTER(create_bar,     tag_detail::fermionic)
    REGISTER(destroy_unbar,  tag_detail::fermionic)
    REGISTER(destroy_bar,    tag_detail::fermionic)
    REGISTER(count_unbar,    tag_detail::bosonic)
    REGISTER(count_bar,      tag_detail::bosonic)
    //REGISTER(e2d,          tag_detail::bosonic)
    //REGISTER(d2e,          tag_detail::bosonic)
    REGISTER(docc,           tag_detail::bosonic)

    #undef REGISTER
    /**********************************************************************/

    // TODO: change term_assistant input paramters
    // For now just pass ident_unbar (will not be used anyway)
    // ATTENTION:  pass fill_unbar --> doesn't affect one-term,
    // but affects 3- and 4-terms --> need to fix it when we use
    // all terms in the hamiltonian!!!
    // NOTE: actually also fill is taken from here and there should be
    // no problem --> need to check again
    chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, dummy, dummy, tag_handler, this->align, this);
    
    std::vector<value_type> & matrix_elements = term_assistant.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);
    
    // TODO: move it up, remove print statement
    int n_pair = lat.size()/2;

    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = term_assistant.idx(m, 0);
        int j = term_assistant.idx(m, 1);
        int k = term_assistant.idx(m, 2);
        int l = term_assistant.idx(m, 3);

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {
            
            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back( boost::make_tuple(0, ident_unbar) );
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                term_descriptor term;
                term.coeff = matrix_elements[m];
                term.push_back( boost::make_tuple(i, count_unbar));
                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.coeff = matrix_elements[m];
                term.push_back( boost::make_tuple(i+n_pair, count_bar));
                this->terms_.push_back(term);
            }
            
            used_elements[m] += 1;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {
            
            // i & j refer to unbarred spinors --> passs create_/destroy/fill_unbar
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill_unbar, matrix_elements[m], i, j, create_unbar, destroy_unbar, tag_handler)
            );
            // i+n_pair, j+n_pair refer to barred spinors --> bar operators
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill_bar, matrix_elements[m], i + n_pair, j + n_pair, create_bar, destroy_bar, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill_unbar, matrix_elements[m], j, i, create_unbar, destroy_unbar, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill_bar, matrix_elements[m], j + n_pair, i + n_pair, create_bar, destroy_bar, tag_handler)
            );
            
            used_elements[m] += 1;
        }

        // Two-electron terms V_ijkl
        // All possible combinations are read in from the FCIDUMP file
        // ijkl stands for unbarred spinors, pqrs stands for barred spinors
        
        // V_iiii
        else if (i == j && j == k && k == l) {
        
        // This case doesn't exist for the relativistic model!
        // We cannot have the same spinor occupied twice!

        /*    
        //------------- 4 unbarred & 0 barred ---------------//
            if (i < n_pair) {
                term_descriptor term;
                term.coeff = matrix_elements[m]*0.5;
                //term.push_back(boost::make_tuple(i, docc));
                term.push_back(boost::make_tuple(i, count_unbar));
                this->terms_.push_back(term);
            }
        //------------- 0 unbarred & 4 barred ---------------//
            else if (i >= n_pair) {
                term_descriptor term;
                term.coeff = matrix_elements[m]*0.5;
                //term.push_back(boost::make_tuple(i, docc));
                term.push_back(boost::make_tuple(i, count_bar));
                this->terms_.push_back(term);
            }
        */    
            used_elements[m] += 1;
        
        }
        
        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            if (i < n_pair && k < n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, i, k, count_unbar, count_unbar);
            } else if (i < n_pair && k >= n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, i, k, count_unbar, count_bar);
            } else if (i >= n_pair && k < n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, i, k, count_bar, count_unbar);
            } else if (i >= n_pair && k >= n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, i, k, count_bar, count_bar);
            }
            
            used_elements[m] += 1;
        }
       
        // V_iikl == V_ijll
        else if ( i==j && j!=k && j!=l || k==l && i!=k && j!=k ) {
            
            int same_idx;
            if (i==j) { same_idx = i; }
            if (k==l) { same_idx = k; }

            // ----------------- ii unbarred ----------------- //
            //
            // ii unbarred, k unbarred, l unbarred
            if (i < n_pair && k < n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii unbarred, k unbarred, l barred
            else if (i < n_pair && k < n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii unbarred, k barred, l unbarred
            else if (i < n_pair && k >= n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii unbarred, kl barred
            else if (i < n_pair && k >= n_pair && l >= n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, same_idx, k, l, create_unbar, destroy_unbar, create_bar, destroy_bar);
                std::cout << "registered: " << i << j << k << l << std::endl;
            }
            // ----------------------------------------------- //
            

            // ------------------ ii barred ------------------ //
            //
            // ii barred, k barred, l barred
            else if (i >= n_pair && k >= n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii barred, k unbarred, l barred
            else if (i >= n_pair && k < n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii barred, k barred, l unbarred
            else if (i >= n_pair && k >= n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ii barred, kl unbarred
            else if (i >= n_pair && k < n_pair && l < n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, same_idx, k, l, create_bar, destroy_bar, create_unbar, destroy_unbar);
                std::cout << "registered: " << i << j << k << l << std::endl;
            }
            // ----------------------------------------------- //
            
            
            // ----------------- ll unbarred ----------------- //
            //
            // i unbarred, j unbarred, ll unbarred
            else if (i < n_pair && j < n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // i unbarred, j barred, ll unbarred
            else if (i < n_pair && j >= n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // i barred, j unbarred, ll unbarred
            else if (i >= n_pair && j < n_pair && l < n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ij barred, ll unbarred
            else if (i >= n_pair && j >= n_pair && l < n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, same_idx, j, i, create_unbar, destroy_unbar, create_bar, destroy_bar);
                std::cout << "registered: " << i << j << k << l << std::endl;
            }
            // ----------------------------------------------- //
            
            
            // ------------------ ll barred ------------------ //
            //
            // i barred, j barred, ll barred
            else if (i >= n_pair && j >= n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // i barred, j unbarred, ll barred
            else if (i >= n_pair && j < n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // i unbarred, j barred, ll barred
            else if (i < n_pair && j >= n_pair && l >= n_pair) {
                // do nothing for the moment
                std::cout << "neglected: " << i << j << k << l << std::endl;
            }
            // ij unbarred, ll barred
            else if (i < n_pair && j < n_pair && l >= n_pair) {
                term_assistant.add_term(this->terms_, matrix_elements[m]*0.5, same_idx, j, i, create_bar, destroy_bar, create_unbar, destroy_unbar);
                std::cout << "registered: " << i << j << k << l << std::endl;
            }
            // ----------------------------------------------- //
            
            used_elements[m] += 1;
        } 

        /*------------- 3 unbarred & 1 barred ---------------*/
        // V_pjkl
        else if (i >= n_pair && j < n_pair && k < n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_unbar, destroy_unbar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // V_iqkl
        else if (i < n_pair && j >= n_pair && k < n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_bar, destroy_unbar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // V_ijrl
        else if (i < n_pair && j < n_pair && k >= n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_unbar, destroy_bar, destroy_unbar);
             
            used_elements[m] += 1;
        }
        
        // V_ijks
        else if (i < n_pair && j < n_pair && k < n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_unbar, destroy_unbar, destroy_bar);
             
            used_elements[m] += 1;
        }

        /*------------- 2 unbarred & 2 barred ---------------*/
        // V_pqkl
        else if (i >= n_pair && j >= n_pair && k < n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_bar, destroy_unbar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // V_pjrl
        else if (i >= n_pair && j < n_pair && k >= n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_unbar, destroy_bar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // V_pjks
        else if (i >= n_pair && j < n_pair && k < n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_unbar, destroy_unbar, destroy_bar);
             
            used_elements[m] += 1;
        }

        // V_iqrl
        else if (i < n_pair && j >= n_pair && k >= n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_bar, destroy_bar, destroy_unbar);
             
            used_elements[m] += 1;
        }


        // V_iqks
        else if (i < n_pair && j >= n_pair && k < n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_bar, destroy_unbar, destroy_bar);
             
            used_elements[m] += 1;
        }

        // V_ijrs
        else if (i < n_pair && j < n_pair && k >= n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_unbar, create_unbar, destroy_bar, destroy_bar);
             
            used_elements[m] += 1;
        }

        /*------------- 1 unbarred & 3 barred ---------------*/
        // V_pqrl
        else if (i >= n_pair && j >= n_pair && k >= n_pair && l < n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_bar, destroy_bar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // V_pqks
        else if (i >= n_pair && j >= n_pair && k < n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_bar, destroy_unbar, destroy_bar);
             
            used_elements[m] += 1;
        }

        // V_pjrs
        else if (i >= n_pair && j < n_pair && k >= n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_unbar, destroy_bar, destroy_bar);
             
            used_elements[m] += 1;
        }
        
        // V_iqrs
        else if (i < n_pair && j >= n_pair && k >= n_pair && l >= n_pair) {
            
            //term_assistant.add_term(this->terms_, i,j,k,l, create_bar, create_bar, destroy_bar, destroy_unbar);
             
            used_elements[m] += 1;
        }

        // --------------------------SEB code------------------------------ //
        /*
        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {

            int same_idx, pos1;

            if      (i==j) { same_idx = i; pos1 = l; }
            else if (k==l) { same_idx = l; pos1 = i; }
            else           { throw std::runtime_error("Term generation logic has failed for V_ijjj term\n"); }

            std::pair<tag_type, value_type> ptag;

            // 1a
            // --> c_l_up * n_i_down * cdag_i_up
            ptag = tag_handler->get_product_tag(count_bar, create_unbar);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, destroy_unbar, tag_handler) );

            // 1a_dagger
            // --> c_i_up * n_i_down * cdag_l_up
            ptag = tag_handler->get_product_tag(destroy_unbar, count_bar);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, -matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, create_unbar, tag_handler) );

            // 1b
            // --> c_l_down * n_i_up * cdag_i_down (1b)
            ptag = tag_handler->get_product_tag(count_unbar, create_bar);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, destroy_bar, tag_handler) );

            // (1b)_dagger
            // --> c_i_down * n_i_up * cdag_l_down
            ptag = tag_handler->get_product_tag(destroy_bar, count_unbar);
            this->terms_.push_back( TermMaker<Matrix, SymmGroup>::positional_two_term(true, fill, -matrix_elements[m] * ptag.second, same_idx, pos1,
                                           ptag.first, create_bar, tag_handler) );

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_unbar, count_unbar);
            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_unbar, count_bar);
            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_bar, count_unbar);
            term_assistant.add_term(this->terms_, matrix_elements[m], i, k, count_bar, count_bar);

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {

            //term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, e2d, d2e);
            //term_assistant.add_term(this->terms_,  matrix_elements[m], i, j, d2e, e2d);
            term_assistant.add_term(this->terms_, -matrix_elements[m], i, j, count_unbar, count_unbar);
            term_assistant.add_term(this->terms_, -matrix_elements[m], i, j, count_bar, count_bar);

            std::pair<tag_type, value_type> ptag1, ptag2;

            // Could insert fill operators without changing the result
            // --> -c_j_up * cdag_j_down * c_i_down * cdag_i_up
            ptag1 = tag_handler->get_product_tag(destroy_bar, create_unbar);
            ptag2 = tag_handler->get_product_tag(destroy_unbar, create_bar);
            term_assistant.add_term(
                this->terms_, -matrix_elements[m] * ptag1.second * ptag2.second, i, j, ptag1.first, ptag2.first
            );

            // --> -c_i_up * cdag_i_down * c_j_down * cdag_j_up
            ptag1 = tag_handler->get_product_tag(destroy_unbar, create_bar);
            ptag2 = tag_handler->get_product_tag(destroy_bar, create_unbar);
            term_assistant.add_term(
                this->terms_, -matrix_elements[m] * ptag1.second * ptag2.second, i, j, ptag1.first, ptag2.first
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
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
            // n_up * cdag_down * c_down <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_unbar, destroy_unbar, create_bar, destroy_bar);
            // n_down * cdag_up * c_up <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_bar, destroy_bar, create_unbar, destroy_unbar);
            // n_down * cdag_down * c_down <--
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, k, l, create_bar, destroy_bar, create_bar, destroy_bar);

            // --> n_up * c_up * cdag_up
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
            // --> n_up * c_down * cdag_down
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_unbar, destroy_unbar, create_bar, destroy_bar);
            // --> n_down * c_up * cdag_up
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_bar, destroy_bar, create_unbar, destroy_unbar);
            // --> n_down * c_down * cdag_down
            term_assistant.add_term(this->terms_, matrix_elements[m], same_idx, l, k, create_bar, destroy_bar, create_bar, destroy_bar);

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

            std::pair<tag_type, value_type> ptag;

            ptag = tag_handler->get_product_tag(create_unbar, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, create_bar , destroy_bar, destroy_unbar
            );
            ptag = tag_handler->get_product_tag(create_bar, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, create_unbar   , destroy_unbar  , destroy_bar
            );
            ptag = tag_handler->get_product_tag(destroy_bar, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, destroy_unbar  , create_unbar   , create_bar
            );
            ptag = tag_handler->get_product_tag(destroy_unbar, fill);
            term_assistant.add_term(
                this->terms_, matrix_elements[m]*ptag.second, same_idx, pos1, pos2, ptag.first, destroy_bar, create_bar , create_unbar
            );

            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_unbar,   destroy_unbar,   create_unbar,   destroy_unbar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_unbar,   destroy_bar, create_bar, destroy_unbar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_bar, destroy_unbar,   create_unbar,   destroy_bar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos1, pos2, create_bar, destroy_bar, create_bar, destroy_bar
            );

            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_unbar,   destroy_unbar,   create_unbar,   destroy_unbar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_unbar,   destroy_bar, create_bar, destroy_unbar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_bar, destroy_unbar,   create_unbar,   destroy_bar
            );
            term_assistant.add_term(
                this->terms_, -matrix_elements[m], same_idx, pos2, pos1, create_bar, destroy_bar, create_bar, destroy_bar
            );

            used_elements[m] += 1;
        }

        // 32 (8x4)-fold degenerate V_ijkl = V_jikl = V_ijlk = V_jilk = V_klij = V_lkij = V_klji = V_lkji * spin
        // V_ijkl -> 24 permutations which fall into 3 equivalence classes of 8 elements (with identical V_ matrix element)
        // coded: 4 index permutations x 4 spin combinations 
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {
            
            // 1
            term_assistant.add_term(this->terms_, i,k,l,j, create_unbar, create_unbar, destroy_unbar, destroy_unbar);
            term_assistant.add_term(this->terms_, i,k,l,j, create_unbar, create_bar, destroy_bar, destroy_unbar);
            term_assistant.add_term(this->terms_, i,k,l,j, create_bar, create_unbar, destroy_unbar, destroy_bar);
            term_assistant.add_term(this->terms_, i,k,l,j, create_bar, create_bar, destroy_bar, destroy_bar);

            // 2
            term_assistant.add_term(this->terms_, i,l,k,j, create_unbar, create_unbar, destroy_unbar, destroy_unbar);
            term_assistant.add_term(this->terms_, i,l,k,j, create_unbar, create_bar, destroy_bar, destroy_unbar);
            term_assistant.add_term(this->terms_, i,l,k,j, create_bar, create_unbar, destroy_unbar, destroy_bar);
            term_assistant.add_term(this->terms_, i,l,k,j, create_bar, create_bar, destroy_bar, destroy_bar);

            // 3
            term_assistant.add_term(this->terms_, j,k,l,i, create_unbar, create_unbar, destroy_unbar, destroy_unbar);
            term_assistant.add_term(this->terms_, j,k,l,i, create_unbar, create_bar, destroy_bar, destroy_unbar);
            term_assistant.add_term(this->terms_, j,k,l,i, create_bar, create_unbar, destroy_unbar, destroy_bar);
            term_assistant.add_term(this->terms_, j,k,l,i, create_bar, create_bar, destroy_bar, destroy_bar);

            // 4
            term_assistant.add_term(this->terms_, j,l,k,i, create_unbar, create_unbar, destroy_unbar, destroy_unbar);
            term_assistant.add_term(this->terms_, j,l,k,i, create_unbar, create_bar, destroy_bar, destroy_unbar);
            term_assistant.add_term(this->terms_, j,l,k,i, create_bar, create_unbar, destroy_unbar, destroy_bar);
            term_assistant.add_term(this->terms_, j,l,k,i, create_bar, create_bar, destroy_bar, destroy_bar);
        
            used_elements[m] += 1;
        }*/
    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    //maquis::cout << bool(it_0 < used_elements.end()) << std::endl;
    //assert( it_0 == used_elements.end() );

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
    
}
    
#endif
