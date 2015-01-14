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

#ifndef REL_QC_MODEL_TS_HPP
#define REL_QC_MODEL_TS_HPP

template <class Matrix, class SymmGroup>
rel_qc_model_ts<Matrix, SymmGroup>::rel_qc_model_ts(Lattice const & lat_, BaseParameters & parms_)
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
					
	// A -> empty, B -> occupied
    typename SymmGroup::charge A(0), B(0);
    B[0]=1;

    // create physical indices
    phys_indices.resize(lat.size());
    for (std::size_t site = 0; site < lat.size(); ++site)
    {
        // Set double group
        int dg = lat.get_prop<int>("irrep", site);
        Index<SymmGroup> loc;
        B[1] = dg;
        loc.insert(std::make_pair(A, 1));
        loc.insert(std::make_pair(B, 1));
        phys_indices[site] = loc;
    }
	
	B[1] = 0;
	
    op_t create_unbar_op, create_bar_op, destroy_unbar_op, destroy_bar_op,
         count_unbar_op, count_bar_op, ident_transfer_op, fill_transfer_op,
         ident_unbar_op, ident_bar_op, fill_unbar_op, fill_bar_op;
	
    ident_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_op.insert_block(Matrix(1, 1, 1), B, B);
	
    create_op.insert_block(Matrix(1, 1, 1), A, B);

    destroy_op.insert_block(Matrix(1, 1, 1), B, A);

    count_op.insert_block(Matrix(1, 1, 1), B, B);

    fill_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_op.insert_block(Matrix(1, 1, -1), B, B);
	
    op_t tmp;
    tag_type dummy;

	// Anticommutation relations!!
 	//gemm(fill_unbar_op, destroy_unbar_op, tmp);
    //destroy_unbar_op = tmp;
 	//gemm(fill_bar_op, destroy_bar_op, tmp);
    //destroy_bar_op = tmp;
	
    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    #define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);

    REGISTER(ident,    tag_detail::bosonic)
    REGISTER(fill,     tag_detail::bosonic)
    REGISTER(create,   tag_detail::fermionic)
    REGISTER(destroy,  tag_detail::fermionic)
    REGISTER(count,    tag_detail::bosonic)

    #undef REGISTER
    /**********************************************************************/

    // TODO: change term_assistant input paramters
    rel_chem_detail::ChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, dummy, dummy, tag_handler, this->align, this);
    
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

		#define onsite
		#ifdef onsite
        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                term_descriptor term;
                term.coeff = matrix_elements[m]*2.0;
                term.push_back( boost::make_tuple(i, count));
                this->terms_.push_back(term);
            }
            {
                term_descriptor term;
                term.coeff = matrix_elements[m]*2.0;
                term.push_back( boost::make_tuple(i+n_pair, count));
                this->terms_.push_back(term);
            }
            
            used_elements[m] += 1;
        }
		#endif

        #define hopping
		#ifdef hopping
        // Hopping term t_ij 
        else if ( i!=j && k == -1 && l == -1) {
            
            //matrix_elements[m] = 0.0;
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, i, j, create, destroy, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, i + n_pair, j + n_pair, create, destroy, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, j, i, create, destroy, tag_handler)
            );
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(this, 
                true, dummy, matrix_elements[m]*2.0, j + n_pair, i + n_pair, create, destroy, tag_handler)
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
        else if ( (i==j && j==k && k!=l) || (i==j && j==l && j!=k) || (i!=j && i==k && k==l) || (i!=j && j==k && k==l) ) {
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
            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if (i < n_pair && k < n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, count, count);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k+n_pair, count, count);
                } else if (i < n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, count, count);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k-n_pair, count, count);
                } else if (i >= n_pair && k < n_pair) {
                } else if (i >= n_pair && k >= n_pair) {
                }
            used_elements[m] += 1;
            }
        }
       
        // V_ijji
        else if ( i==l && j==k && i!=j) {
            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if (i < n_pair && j < n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, j, count, count);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, j+n_pair, count, count);
                } else if (i < n_pair && j >= n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, j, count, count);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, j-n_pair, count, count);
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
            
            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if        (i < n_pair && k <  n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, l, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k+n_pair, l+n_pair, create, destroy, create, destroy);
                } else if (i < n_pair && k >= n_pair && l >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i, k, l, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											i+n_pair, k-n_pair, l-n_pair, create, destroy, create, destroy);
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

            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && k  < n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k, i, j, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k+n_pair, i+n_pair, j+n_pair, create, destroy, create, destroy);
                } else if (i <  n_pair && j <  n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k, i, j, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, matrix_elements[m],
											k-n_pair, i+n_pair, j+n_pair, create, destroy, create, destroy);
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

            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && k <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, k, j, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, k+n_pair, j+n_pair, create, destroy, create, destroy);
                } else if (i <  n_pair && j >= n_pair && k >= n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i, k, j, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											i+n_pair, k-n_pair, j-n_pair, create, destroy, create, destroy);
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

            //matrix_elements[m] = 0.0;
            if (is_term_allowed(i,j,k,l)) {
                if        (i <  n_pair && j <  n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j, i, l, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j+n_pair, i+n_pair, l+n_pair, create, destroy, create, destroy);
                } else if (i <  n_pair && j >= n_pair && l <  n_pair) {
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j, i, l, create, destroy, create, destroy);
                    term_assistant.add_term(this->terms_, -matrix_elements[m],
											j-n_pair, i+n_pair, l+n_pair, create, destroy, create, destroy);
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
        else if (i!=j && j!=k && k!=l && i!=k && j!=l && i!=l) {
            
            if (is_term_allowed(i,j,k,l)) {
// 				maquis::cout << "4_term" << std::endl;
                if        (i <  n_pair && j <  n_pair && k <  n_pair && l <  n_pair) {

					//matrix_elements[m] = 0.0;
					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
											create, create, destroy, destroy);
											
                    term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l+n_pair,j+n_pair,
											create, create, destroy, destroy);
                
				} else if (i <  n_pair && j <  n_pair && k >= n_pair && l >= n_pair) {
                
                    //matrix_elements[m] = 0.0;
					//matrix_elements[m] *= 10000.0;
 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
 											create, create, destroy, destroy);
											
                    term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l-n_pair,j+n_pair,
 											create, create, destroy, destroy);
											
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
					
					//maquis::cout << "##############################################\n";
					//maquis::cout << "(" << i << "," << j << "," << k << "," << l << ")" << std::endl;
					//maquis::cout << "(" << i+n_pair << "," << j-n_pair<< "," << k+n_pair << "," << l-n_pair << ")" << std::endl;
					//if ((i == 2 && j == 5 && k == 0 && l == 6)){ //|| (i == 0 && j == 5 && k == 2 && l == 6) || (i == 0 && j == 6 && k == 2 && l == 5) || (i == 2 && j == 6 && k == 0 && l == 5) ) {
						//maquis::cout << "(" << i << "," << j << "," << k << "," << l << "): " << 2.0*matrix_elements[m] << std::endl;
					//matrix_elements[m] *= 10000.0;
 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
 											create, create, destroy, destroy);
					//}
					//if ((i == 2 && j == 4 && k == 1 && l == 6)) { //|| (i == 1 && j == 4 && k == 2 && l == 6) || (i == 1 && j == 6 && k == 2 && l == 4) || (i == 2 && j == 6 && k == 1 && l == 4)) {
					//	//maquis::cout << "(" << i+n_pair << "," << j-n_pair<< "," << k+n_pair << "," << l-n_pair << "): " << 2.0*matrix_elements[m] << std::endl;	
					//matrix_elements[m] *= 10000.0;
 					term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k+n_pair,l-n_pair,j-n_pair,
 											create, create, destroy, destroy);
					//}
                } else if (i <  n_pair && j >= n_pair && k >= n_pair && l <  n_pair) {
					
                    //matrix_elements[m] *= 0.0;
                    term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i,k,l,j,
 											create, create, destroy, destroy);
											
                    term_assistant.add_term(this->terms_, matrix_elements[m],n_pair,i+n_pair,k-n_pair,l+n_pair,j-n_pair,
 											create, create, destroy, destroy);
											
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
    //assert( it_0 == used_elements.end() );

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}
    
#endif
