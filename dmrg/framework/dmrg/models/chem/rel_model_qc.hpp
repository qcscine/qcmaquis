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
rel_qc_model<Matrix, SymmGroup>::rel_qc_model(Lattice const & lat_, BaseParameters & parms_)
: lat(lat_)
, parms(parms_)
, tag_handler(new table_type())
{
    // find the highest irreducible representation number 
    // used to generate ops for all irreps 0..max_irrep
    // TODO: obsolete, probably not used anywhere -> remove it
    max_irrep = 0;
    for (pos_t p=0; p < lat.size(); ++p)
        max_irrep = (lat.get_prop<typename SymmGroup::subcharge>("irrep", p) > max_irrep)
                        ? lat.get_prop<typename SymmGroup::subcharge>("irrep", p) : max_irrep;

	// Hardcoded for the time being
	max_irrep = 128;
					
	// A -> empty, B -> unbarred spinor, C-> barred spinor
    typename SymmGroup::charge A(0), B(0), C(0);
    B[0]=1; C[1]=1;

    // create physical indices
    phys_indices.resize(lat.size());
    for (std::size_t site = 0; site < lat.size(); ++site)
    {
        // Set double group
        int dg = lat.get_prop<int>("irrep", site);

        Index<SymmGroup> loc;
        if (site < lat.size()/2)
        {   // create unbarred bases in the first half
			B[2] = dg;
            loc.insert(std::make_pair(A, 1));
            loc.insert(std::make_pair(B, 1));
        }
        else
        {   // barred bases in the second half
			C[2] = dg;
            loc.insert(std::make_pair(A, 1));
            loc.insert(std::make_pair(C, 1));
        }
        phys_indices[site] = loc;
    }
	
	B[2] = 0;
	C[2] = 0;
	
    op_t create_unbar_op, create_bar_op, destroy_unbar_op, destroy_bar_op,
         count_unbar_op, count_bar_op,
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

    fill_unbar_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_unbar_op.insert_block(Matrix(1, 1, -1), B, B);
	
    fill_bar_op.insert_block(Matrix(1, 1, 1), A, A);
	fill_bar_op.insert_block(Matrix(1, 1, -1), C, C);

    op_t tmp;
    tag_type dummy;

	// Anticommutation relations!!
// 	gemm(create_bar_op, fill_bar_op, tmp);
//     create_bar_op = tmp;
// 	gemm(fill_bar_op, destroy_bar_op, tmp);
//     destroy_bar_op = tmp;

// 	maquis::cout << create_unbar_op << std::endl;
// 	maquis::cout << fill_unbar_op << std::endl;
	
// 	gemm(create_unbar_op, fill_unbar_op, tmp);
	
// 	create_unbar_op = tmp;
	
// 	maquis::cout << create_unbar_op << std::endl;
	
// 	gemm(fill_unbar_op, destroy_unbar_op, tmp);
// 	destroy_unbar_op = tmp;
	
// 	maquis::cout << create_unbar_op << std::endl;
// 	maquis::cout << create_bar_op << std::endl;
// 	maquis::cout << destroy_unbar_op << std::endl;
// 	maquis::cout << destroy_bar_op << std::endl;
	
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

    #undef REGISTER
    /**********************************************************************/

	//#define check_op
	#ifdef check_op
	for (int ii=0; ii < 10; ++ii) {
		op_t debug_op = tag_handler->get_op(ii);
		maquis::cout << debug_op << std::endl;
	}
	#endif
	
	// DEBUG_TAGHandler
// 	#define DEBUG_TAGHandler
	#ifdef DEBUG_TAGHandler
	maquis::cout << "DEBUG_TAGHandler" << std::endl;
	maquis::cout << "fill_unbar: " << fill_unbar << std::endl << tag_handler->get_op(fill_unbar) << std::endl;
	maquis::cout << "destroy_unbar: " << destroy_unbar << std::endl << tag_handler->get_op(destroy_unbar) << std::endl;
	std::pair<tag_type,value_type> debug_prod_op = tag_handler->get_product_tag(fill_unbar,destroy_unbar);
	maquis::cout << "get_product_tag(fill_unbar,destroy_unbar): " << debug_prod_op.first << std::endl;
	maquis::cout << tag_handler->get_op(debug_prod_op.first) << std::endl;
	maquis::cout << debug_prod_op.second << std::endl;
	#endif
	
    // TODO: change term_assistant input paramters
    chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, dummy, dummy, tag_handler, this->align, this);
    
    std::vector<value_type> & matrix_elements = term_assistant.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);
    
    // TODO: move it up
    int n_pair = lat.size()/2;

    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = term_assistant.idx(m, 0);
        int j = term_assistant.idx(m, 1);
        int k = term_assistant.idx(m, 2);
        int l = term_assistant.idx(m, 3);

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {
            
            term_descriptor term;
            term.coeff = matrix_elements[m]*2.0;
            term.push_back( boost::make_tuple(0, ident_unbar) );
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                term_descriptor term;
                term.coeff = matrix_elements[m]*2.0;
                term.push_back( boost::make_tuple(i, count_unbar));
                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.coeff = matrix_elements[m]*2.0;
                term.push_back( boost::make_tuple(i+n_pair, count_bar));
                this->terms_.push_back(term);
            }
            
            used_elements[m] += 1;
        }

		#define hopping
		#ifdef hopping
        // Hopping term t_ij 
        else if (k == -1 && l == -1) {
            
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, i, j, create_unbar, destroy_unbar, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, i + n_pair, j + n_pair, create_bar, destroy_bar, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, j, i, create_unbar, destroy_unbar, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, j + n_pair, i + n_pair, create_bar, destroy_bar, tag_handler)
            );
            
            used_elements[m] += 1;
        }
        #endif

        // V_iiii
        // This case doesn't exist for the relativistic model!
        else if (i == j && j == k && k == l) {
            used_elements[m] += 1;
        }
        
        // V_ijjj & V_jijj & V_jjij & V_jjji
        // This case doesn't exist for the relativistic model!
        else if ( (i==j && j==k && k!=l) || (i==j && j!=k && j==l) || (i!=j && i==k && k==l) || (i!=j && j==k && k==l) ) {
            used_elements[m] += 1;
        }
        
        // V_ijij
        // This case doesn't exist for the relativistic model!
        else if ( i==k && j==l && i!=j) {
            used_elements[m] += 1;
        }

		#define all_2_terms
		#ifdef all_2_terms
        // V_iijj
        else if ( i==j && k==l && j!=k) {
            
            if (is_term_allowed(i,j,k,l)) {
                if (i < n_pair && k < n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, count_unbar, count_unbar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k+n_pair, count_bar, count_bar);
                } else if (i < n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, count_unbar, count_bar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k-n_pair, count_bar, count_unbar);
                } else if (i >= n_pair && k < n_pair) {
                } else if (i >= n_pair && k >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }
       
        // V_ijji
        else if ( i==l && j==k && i!=j) {
            
            if (is_term_allowed(i,j,k,l)) {
                if (i < n_pair && j < n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, j, count_unbar, count_unbar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, j+n_pair, count_bar, count_bar);
                } else if (i < n_pair && j >= n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, j, count_unbar, count_bar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, j-n_pair, count_bar, count_unbar);
                } else if (i >= n_pair && j < n_pair) {
                } else if (i >= n_pair && j >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }
        
        #endif
		
		#define all_3_terms
		#ifdef all_3_terms
        // V_iikl
        else if ( i==j && j!=k && j!=l ) {
            
            if (is_term_allowed(i,j,k,l)) {
                if        (i < n_pair && k <  n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, l, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k+n_pair, l+n_pair, create_bar, destroy_bar, create_bar, destroy_bar);
                } else if (i < n_pair && k >= n_pair && l >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, l, create_unbar, destroy_unbar, create_bar, destroy_bar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k-n_pair, l-n_pair, create_bar, destroy_bar, create_unbar, destroy_unbar);
                } else if (i < n_pair && k >= n_pair && l <  n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											i, k, l, create_unbar, destroy_unbar, create_bar, destroy_unbar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											i+n_pair, k-n_pair, l+n_pair, create_bar, destroy_bar, create_unbar, destroy_bar);
                } else if (i < n_pair && k <  n_pair && l >= n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											i, k, l, create_unbar, destroy_unbar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											i+n_pair, k+n_pair, l-n_pair, create_bar, destroy_bar, create_bar, destroy_unbar);
                } else if (i >= n_pair && k <  n_pair && l <  n_pair) {
                } else if (i >= n_pair && k >= n_pair && l <  n_pair) {
                } else if (i >= n_pair && k <  n_pair && l >= n_pair) {
                } else if (i >= n_pair && k >= n_pair && l >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }
         
        // V_ijkk   
        else if (k==l && i!=k && j!=k ) {

            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && k  < n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k, i, j, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k+n_pair, i+n_pair, j+n_pair, create_bar, destroy_bar, create_bar, destroy_bar);
                } else if (i <  n_pair && j <  n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k, i, j, create_bar, destroy_bar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k-n_pair, i+n_pair, j+n_pair, create_unbar, destroy_unbar, create_bar, destroy_bar);
                } else if (i <  n_pair && j >= n_pair && k <  n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											k, i, j, create_unbar, destroy_unbar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											k+n_pair, i+n_pair, j-n_pair, create_bar, destroy_bar, create_bar, destroy_unbar);
                } else if (i <  n_pair && j >= n_pair && k >= n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											k, i, j, create_bar, destroy_bar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],
// 											k-n_pair, i+n_pair, j-n_pair, create_unbar, destroy_unbar, create_bar, destroy_unbar);
                } else if (i >= n_pair && j >= n_pair && k >= n_pair) {
                } else if (i >= n_pair && j >= n_pair && k <  n_pair) {
                } else if (i >= n_pair && j <  n_pair && k >= n_pair) {
                } else if (i >= n_pair && j <  n_pair && k <  n_pair) {
                }
            used_elements[m] += 1;
            }
        }

        // V_ijil & V_ijkj
        // This case doesn't exist for the relativistic model!
        else if ( i==k && i!=j && i!=l || j==l && j!=i && j!=k ) {
            used_elements[m] += 1;
        }
        
        // V_ijki
        else if ( i!=j && j!=k && i==l ) {

            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && k <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, k, j, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, k+n_pair, j+n_pair, create_bar, destroy_bar, create_bar, destroy_bar);
                } else if (i <  n_pair && j >= n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, k, j, create_unbar, destroy_unbar, create_bar, destroy_bar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, k-n_pair, j-n_pair, create_bar, destroy_bar, create_unbar, destroy_unbar);
                } else if (i <  n_pair && j >= n_pair && k <  n_pair) {
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											i, k, j, create_unbar, destroy_unbar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											i+n_pair, k+n_pair, j-n_pair, create_bar, destroy_bar, create_bar, destroy_unbar);
                } else if (i <  n_pair && j <  n_pair && k >= n_pair) {
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											i, k, j, create_unbar, destroy_unbar, create_bar, destroy_unbar);
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											i+n_pair, k-n_pair, j+n_pair, create_bar, destroy_bar, create_unbar, destroy_bar);
                } else if (i >= n_pair && j <  n_pair && k <  n_pair) {
                } else if (i >= n_pair && j >= n_pair && k <  n_pair) {
                } else if (i >= n_pair && j <  n_pair && k >= n_pair) {
                } else if (i >= n_pair && j >= n_pair && k >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }

        // V_ijjl
        else if ( i!=j && j==k && k!=l ) {

            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j, i, l, create_unbar, destroy_unbar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j+n_pair, i+n_pair, l+n_pair, create_bar, destroy_bar, create_bar, destroy_bar);
                } else if (i <  n_pair && j >= n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j, i, l, create_bar, destroy_bar, create_unbar, destroy_unbar);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j-n_pair, i+n_pair, l+n_pair, create_unbar, destroy_unbar, create_bar, destroy_bar);
                } else if (i <  n_pair && j <  n_pair && l >=  n_pair) {
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											j, i, l, create_unbar, destroy_unbar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											j+n_pair, i+n_pair, l-n_pair, create_bar, destroy_bar, create_bar, destroy_unbar);
                } else if (i <  n_pair && j >= n_pair && l >= n_pair) {
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											j, i, l, create_bar, destroy_bar, create_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, -matrix_elements[m],
// 											j-n_pair, i+n_pair, l-n_pair, create_unbar, destroy_unbar, create_bar, destroy_unbar);
                } else if (i >= n_pair && j <  n_pair && l <  n_pair) {
                } else if (i >= n_pair && j >= n_pair && l <  n_pair) {
                } else if (i >= n_pair && j <  n_pair && l >= n_pair) {
                } else if (i >= n_pair && j >= n_pair && l >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }
		#endif
		
		#define all_4_terms
		#ifdef all_4_terms
        // V_ijkl
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {
            
            if (is_term_allowed(i,j,k,l)) {
// 				maquis::cout << "4_term" << std::endl;
                if        (i <  n_pair && j <  n_pair && k <  n_pair && l <  n_pair) {
                    
					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
											create_unbar, create_unbar, destroy_unbar, destroy_unbar);
											
                    term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l+n_pair,j+n_pair,
											create_bar, create_bar, destroy_bar, destroy_bar);
                
				} else if (i <  n_pair && j <  n_pair && k >= n_pair && l >= n_pair) {
                
// 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_bar, destroy_bar, destroy_unbar);
											
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l-n_pair,j+n_pair,
// 											create_bar, create_unbar, destroy_unbar, destroy_bar);
											
                } else if (i <  n_pair && j <  n_pair && k >= n_pair && l <  n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_bar, destroy_unbar, destroy_unbar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l+n_pair,j+n_pair,
// 											create_bar, create_unbar, destroy_bar, destroy_bar);
                } else if (i <  n_pair && j <  n_pair && k <  n_pair && l >= n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_unbar, destroy_bar, destroy_unbar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l-n_pair,j+n_pair,
// 											create_bar, create_bar, destroy_unbar, destroy_bar);
                } else if (i <  n_pair && j >= n_pair && k <  n_pair && l >= n_pair) {
					
// 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_unbar, destroy_bar, destroy_bar);
											
// 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l-n_pair,j-n_pair,
// 											create_bar, create_bar, destroy_unbar, destroy_unbar);
											
                } else if (i <  n_pair && j >= n_pair && k >= n_pair && l <  n_pair) {
					
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_bar, destroy_unbar, destroy_bar);
											
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l+n_pair,j-n_pair,
// 											create_bar, create_unbar, destroy_bar, destroy_unbar);
											
                } else if (i <  n_pair && j >= n_pair && k <  n_pair && l <  n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_unbar, destroy_unbar, destroy_bar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l+n_pair,j-n_pair,
// 											create_bar, create_bar, destroy_bar, destroy_unbar);
                } else if (i <  n_pair && j >= n_pair && k >= n_pair && l >= n_pair) {
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
// 											create_unbar, create_bar, destroy_bar, destroy_bar);
//                     term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l-n_pair,j-n_pair,
// 											create_bar, create_unbar, destroy_unbar, destroy_unbar);
                } else if (i >= n_pair && j >= n_pair && k >= n_pair && l >= n_pair) {
                } else if (i >= n_pair && j >= n_pair && k <  n_pair && l <  n_pair) {
                } else if (i >= n_pair && j >= n_pair && k <  n_pair && l >= n_pair) {
                } else if (i >= n_pair && j >= n_pair && k >= n_pair && l <  n_pair) {
                } else if (i >= n_pair && j <  n_pair && k >= n_pair && l <  n_pair) {
                } else if (i >= n_pair && j <  n_pair && k <  n_pair && l >= n_pair) {
                } else if (i >= n_pair && j <  n_pair && k >= n_pair && l >= n_pair) {
                } else if (i >= n_pair && j <  n_pair && k <  n_pair && l <  n_pair) {
                }
            used_elements[m] += 1;
            }
        }
        #endif

    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    maquis::cout << bool(it_0 < used_elements.end()) << std::endl;
//     assert( it_0 == used_elements.end() );

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}
    
#endif
