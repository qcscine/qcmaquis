/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef REL_QC_MODEL_HPP
#define REL_QC_MODEL_HPP

template <class SymmGroup>
rel_qc_model<SymmGroup>::rel_qc_model(Lattice const & lat_, BaseParameters & parms_)
    : lat(lat_), parms(parms_), tag_handler(new table_type())
{
    // Initialize double group table
    SymmGroup::initialize_dg_table(parms);

	// A -> empty, B -> occupied
    typename SymmGroup::charge A(0), B(0);
    B[0]=1;

    // create physical indices
    //phys_indices.resize(SymmGroup::get_max_irrep());
    for (std::size_t irrep = 0; irrep < SymmGroup::get_max_irrep()+1; ++irrep)
    {
        Index<SymmGroup> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(B, irrep),1));

        phys_indices.push_back(phys);
    }

    op_t create_op, destroy_op, count_op,
         ident_op, fill_op;

    ident_op.insert_block(Matrix(1, 1, 1), A, A);
    ident_op.insert_block(Matrix(1, 1, 1), B, B);

    create_op.insert_block(Matrix(1, 1, 1), A, B);

    destroy_op.insert_block(Matrix(1, 1, 1), B, A);

    count_op.insert_block(Matrix(1, 1, 1), B, B);

    fill_op.insert_block(Matrix(1, 1, 1), A, A);
    fill_op.insert_block(Matrix(1, 1, -1), B, B);

    #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

    GENERATE_SITE_SPECIFIC(ident_op)
    GENERATE_SITE_SPECIFIC(fill_op)
    GENERATE_SITE_SPECIFIC(create_op)
    GENERATE_SITE_SPECIFIC(destroy_op)
    GENERATE_SITE_SPECIFIC(count_op)

    #undef GENERATE_SITE_SPECIFIC

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    #define REGISTER(op, kind) op = this->register_site_specific(op ## _ops, kind);

    REGISTER(ident,    tag_detail::bosonic)
    REGISTER(fill,     tag_detail::bosonic)
    REGISTER(create,   tag_detail::fermionic)
    REGISTER(destroy,  tag_detail::fermionic)
    REGISTER(count,    tag_detail::bosonic)

    #undef REGISTER

    /**********************************************************************/
}

template<class SymmGroup>
void rel_qc_model<SymmGroup>::create_terms()
{
    chem::detail::RelChemHelper<Matrix, SymmGroup> term_assistant(parms, lat, ident, fill, tag_handler);
    auto& matrix_elements = term_assistant.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);

    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = term_assistant.idx(m, 0);
        int j = term_assistant.idx(m, 1);
        int k = term_assistant.idx(m, 2);
        int l = term_assistant.idx(m, 3);
        
        auto matrixElement = static_cast<value_type>(matrix_elements[m]);

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {

            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back( std::make_pair(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]) );
            this->terms_.push_back(term);
            used_elements[m] += 1;
        }

		#define onsite
		#ifdef onsite
        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {
            {
                term_descriptor term;
                term.coeff = matrixElement;
                term.push_back( std::make_pair(i, count[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
                this->terms_.push_back(term);
            }
            used_elements[m] += 1;
        }
		#endif

        #define hopping
		#ifdef hopping
        // Hopping term t_ij
        else if ( i!=j && k == -1 && l == -1) {
            this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                true, fill, matrixElement, i, j, create, destroy, tag_handler, lat)
            );
                if (!parms["MAGNETIC"]) { // with an external magnetic field applied - Kramers symmetry is broken
                  this->terms_.push_back(TermMaker<Matrix, SymmGroup>::positional_two_term(
                      true, fill, chem::detail::cconj(matrix_elements[m]), j, i, create, destroy, tag_handler, lat)
                  );
                }
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
            if (is_term_allowed(i,j,k,l)) {
                term_assistant.add_term(this->terms_, matrixElement, i, k, count, count);
            	used_elements[m] += 1;
            }
        }

        // V_ijji
        else if ( i==l && j==k && i!=j) {
            if (is_term_allowed(i,j,k,l)) {
                term_assistant.add_term(this->terms_, -matrixElement, i, j, count, count);
            	used_elements[m] += 1;
            }
        }
        #endif

 		#define all_3_terms
		#ifdef all_3_terms
        // V_iikl
        else if ( i==j && j!=k && j!=l ) {
            if (is_term_allowed(i,j,k,l)) {
                term_assistant.add_term(this->terms_, matrixElement, i, k, l, create, destroy, create, destroy);
            	used_elements[m] += 1;
            }
        }

        // V_ijkk
        else if (k==l && i!=k && j!=k ) {
            if (is_term_allowed(i,j,k,l)) {
                term_assistant.add_term(this->terms_, matrixElement, k, i, j, create, destroy, create, destroy);
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
                term_assistant.add_term(this->terms_, -matrixElement, i, k, j, create, destroy, create, destroy);
            	used_elements[m] += 1;
            }
        }

        // V_ijjl
        else if ( i!=j && j==k && k!=l ) {
            if (is_term_allowed(i,j,k,l)) {
                term_assistant.add_term(this->terms_, -matrixElement, j, i, l, create, destroy, create, destroy);
            	used_elements[m] += 1;
            }
        }
		#endif

 		#define all_4_terms
		#ifdef all_4_terms
        // V_ijkl
        else if (i!=j && j!=k && k!=l && i!=k && j!=l && i!=l) {
            if (is_term_allowed(i,j,k,l)) {
				term_assistant.add_term(this->terms_, matrixElement, i, k, l, j, create, create, destroy, destroy);
            	used_elements[m] += 1;
            }
        }
        #endif

    } // matrix_elements for

    // make sure all elements have been used
    std::vector<int>::iterator it_0;
    it_0 = std::find(used_elements.begin(), used_elements.end(), 0);
    assert( it_0 == used_elements.end() );

    term_assistant.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
