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

#ifndef QC_SU2_HPP
#define QC_SU2_HPP


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
    create_fill_op.spin() = one_half_up;
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_fill_op.insert_block(Matrix(1,1,1), D, B);
    create_fill_op.insert_block(Matrix(1,1,1), D, C);

    op_t destroy_op;
    destroy_op.spin() = one_half_down;
    destroy_op.insert_block(Matrix(1,1,1), A, B);
    destroy_op.insert_block(Matrix(1,1,1), A, C);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
    destroy_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

    op_t destroy_fill_op;
    destroy_fill_op.spin() = one_half_up;
    destroy_fill_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_op.insert_block(Matrix(1,1,1), A, C);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), B, D);
    destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

    op_t create_op;
    create_op.spin() = one_half_down;
    create_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
    create_op.insert_block(Matrix(1,1,-1), D, B);
    create_op.insert_block(Matrix(1,1,-1), D, C);

    /*************************************************************/

    op_t create_fill_couple_down_op = create_fill_op;
    create_fill_couple_down_op.spin() = one_half_down;

    op_t destroy_fill_couple_down_op = destroy_fill_op;
    destroy_fill_couple_down_op.spin() = one_half_down;

    op_t create_couple_up_op = create_op;
    create_couple_up_op.spin() = one_half_up;

    op_t destroy_couple_up_op = destroy_op;
    destroy_couple_up_op.spin() = one_half_up;

    /*************************************************************/

    op_t create_fill_count_op;
    create_fill_count_op.spin() = one_half_up;
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
    create_fill_count_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);

    op_t destroy_count_op;
    destroy_count_op.spin() = one_half_down;
    destroy_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t destroy_fill_count_op;
    destroy_fill_count_op.spin() = one_half_up;
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, B);
    destroy_fill_count_op.insert_block(Matrix(1,1,1), A, C);

    op_t create_count_op;
    create_count_op.spin() = one_half_down;
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
    flip_to_S2_op.spin() = one_up;
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2)), B, B);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)), C, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  B, C);
    flip_to_S2_op.insert_block(Matrix(1,1,std::sqrt(3./2.)),  C, B);

    op_t flip_to_S0_op = flip_to_S2_op;
    flip_to_S0_op.spin() = one_down;

    op_t flip_S0_op = flip_to_S2_op;
    flip_S0_op.spin() = one_flat;

    /**********************************************************************/
    /*** Create operator tag table ****************************************/
    /**********************************************************************/

    #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = this->generate_site_specific_ops(opname);

    GENERATE_SITE_SPECIFIC(ident_op)
    GENERATE_SITE_SPECIFIC(ident_full_op)
    GENERATE_SITE_SPECIFIC(fill_op)

    GENERATE_SITE_SPECIFIC(create_fill_op)
    GENERATE_SITE_SPECIFIC(create_op)
    GENERATE_SITE_SPECIFIC(destroy_fill_op)
    GENERATE_SITE_SPECIFIC(destroy_op)

    GENERATE_SITE_SPECIFIC(create_fill_couple_down_op)
    GENERATE_SITE_SPECIFIC(destroy_fill_couple_down_op)
    GENERATE_SITE_SPECIFIC(create_couple_up_op)
    GENERATE_SITE_SPECIFIC(destroy_couple_up_op)

    GENERATE_SITE_SPECIFIC(create_fill_count_op)
    GENERATE_SITE_SPECIFIC(create_count_op)
    GENERATE_SITE_SPECIFIC(destroy_fill_count_op)
    GENERATE_SITE_SPECIFIC(destroy_count_op)

    GENERATE_SITE_SPECIFIC(count_op)
    GENERATE_SITE_SPECIFIC(docc_op)
    GENERATE_SITE_SPECIFIC(e2d_op)
    GENERATE_SITE_SPECIFIC(d2e_op)
    GENERATE_SITE_SPECIFIC(flip_S0_op)
    GENERATE_SITE_SPECIFIC(flip_to_S2_op)
    GENERATE_SITE_SPECIFIC(flip_to_S0_op)
    GENERATE_SITE_SPECIFIC(count_fill_op)

    #define REGISTER(op, kind) op = this->register_site_specific(op ## _ops, kind);

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
}

template <class Matrix, class SymmGroup>
void qc_su2<Matrix, SymmGroup>::create_terms()
{
    typedef typename SymmGroup::subcharge subcharge;
    subcharge N = SymmGroup::particleNumber(this->total_quantum_numbers(parms));

    /*************************************************************/
    typename TermMakerSU2<Matrix, SymmGroup>::OperatorBundle create_pkg, destroy_pkg;
    typename TermMakerSU2<Matrix, SymmGroup>::OperatorBundle create_count_pkg, destroy_count_pkg;

    create_pkg.couple_up = create_couple_up;
    create_pkg.couple_down = create;
    create_pkg.fill_couple_up = create_fill;
    create_pkg.fill_couple_down = create_fill_couple_down;

    destroy_pkg.couple_up = destroy_couple_up;
    destroy_pkg.couple_down = destroy;
    destroy_pkg.fill_couple_up = destroy_fill;
    destroy_pkg.fill_couple_down = destroy_fill_couple_down;

    create_count_pkg.couple_down = create_count;
    create_count_pkg.fill_couple_up = create_fill_count;

    destroy_count_pkg.couple_down = destroy_count;
    destroy_count_pkg.fill_couple_up = destroy_fill_count;

    /**********************************************************************/

    op_collection.ident     .no_couple = ident;
    op_collection.ident_full.no_couple = ident_full;
    op_collection.fill      .no_couple = fill;

    op_collection.create               = create_pkg;
    op_collection.destroy              = destroy_pkg;

    op_collection.count     .no_couple = count;
    op_collection.count     .fill_no_couple = count_fill;

    op_collection.create_count         = create_count_pkg;       
    op_collection.destroy_count        = destroy_count_pkg;       

    op_collection.e2d       .no_couple = e2d;
    op_collection.d2e       .no_couple = d2e;
    op_collection.docc      .no_couple = docc;

    op_collection.flip      .no_couple = flip_S0;
    op_collection.flip      .couple_up = flip_to_S2;
    op_collection.flip      .couple_down = flip_to_S0;

    /**********************************************************************/

    chem_detail::ChemHelperSU2<Matrix, SymmGroup> ta(parms, lat, tag_handler);
    alps::numeric::matrix<Lattice::pos_t> idx_ = ta.getIdx();
    std::vector<value_type> matrix_elements = ta.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);

    typedef TermMakerSU2<Matrix, SymmGroup> TM;

    /**********************************************************************/

    using boost::lambda::_1;
    using boost::bind;
    using chem_detail::ChemHelperSU2;
    using chem_detail::append;
 
    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = idx_(m, 0);
        int j = idx_(m, 1);
        int k = idx_(m, 2);
        int l = idx_(m, 3);

        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back( boost::make_tuple(0, ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]) );
            this->terms_.push_back(term);
            
            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back( boost::make_tuple(i, count[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
            this->terms_.push_back(term);

            used_elements[m] += 1;
            continue;
        }

        // Hopping term t_ij 
        else if (k == -1 && l == -1) {

            // one-electron problems need special attention
            if (N == 1) {

            // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
                this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                    true, ident, std::sqrt(2.)*matrix_elements[m],i,j,create, create_fill, destroy, destroy_fill, lat
                ));
                this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                    true, ident, std::sqrt(2.)*matrix_elements[m],j,i,create, create_fill, destroy, destroy_fill, lat
                ));
            }

            else {

                typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
                typedef std::vector<term_descriptor> term_vec;

                term_vec & vec = this->terms_;

                term_vec terms;
                for (pos_t kk = 0; kk < lat.size(); ++kk)
                {   
                    if (kk == j || kk == i) continue;
                    append(terms, SSUM::three_term(matrix_elements[m] * (1./(N-1)), i,kk,kk,j, op_collection, lat));
                    append(terms, SSUM::three_term(matrix_elements[m] * (1./(N-1)), j,kk,kk,i, op_collection, lat));
                }

                std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_3term, &ta, vec, _1));

                terms.clear();

                append(terms, SSUM::V_term(matrix_elements[m] * (1./(N-1)), i,i,i,j, op_collection, lat));
                append(terms, SSUM::V_term(matrix_elements[m] * (1./(N-1)), j,i,i,i, op_collection, lat));

                append(terms, SSUM::V_term(matrix_elements[m] * (1./(N-1)), i,j,j,j, op_collection, lat));
                append(terms, SSUM::V_term(matrix_elements[m] * (1./(N-1)), j,j,j,i, op_collection, lat));
                std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_2term, &ta, vec, _1));
            }

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {

            term_descriptor term;
            term.coeff = matrix_elements[m];
            term.push_back(boost::make_tuple(i, docc[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // V_ijjj = V_jijj = V_jjij = V_jjji
        else if ( (i==j && j==k && k!=l) || (i!=j && j==k && k==l) ) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            int s, p;
            if      (i==j) { s = i; p = l; }
            else if (k==l) { s = l; p = i; }
            else           { throw std::runtime_error("Term generation logic has failed for V_ijjj term\n"); }

            term_vec & vec = this->terms_;

            term_vec terms;
            append(terms, SSUM::two_term(matrix_elements[m], s,s,s,p, op_collection, lat));
            append(terms, SSUM::two_term(matrix_elements[m], s,p,s,s, op_collection, lat));
            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_2term, &ta, vec, _1));

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            term_vec terms = SSUM::two_term(matrix_elements[m], i,k,k,i, op_collection, lat);
            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_2term, &ta, vec, _1));

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {
            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            term_vec terms;
            append(terms, SSUM::two_term(0.5*matrix_elements[m], i,i,j,j, op_collection, lat));
            append(terms, SSUM::two_term(0.5*matrix_elements[m], j,j,i,i, op_collection, lat));

            append(terms, SSUM::two_term(matrix_elements[m], i,j,i,j, op_collection, lat));

            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_2term, &ta, vec, _1));

            used_elements[m] += 1;
        }

        // 9987 9877

        // 8 (4x2)-fold degenerate V_iilk == V_iikl = V_lkii = V_klii  <--- coded
        //                         V_ijkk == V_jikk = V_kkij = V_kkji  <--- contained above
        else if ( (i==j && j!=k && k!=l) || (k==l && i!=j && j!=k)) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            term_vec terms;
            if (i==j)
            {
                append(terms, SSUM::three_term(matrix_elements[m], i,k,l,i, op_collection, lat));
                append(terms, SSUM::three_term(matrix_elements[m], i,l,k,i, op_collection, lat));
            }
            else // (k==l)
            {
                append(terms, SSUM::three_term(matrix_elements[m], i,k,k,j, op_collection, lat));
                append(terms, SSUM::three_term(matrix_elements[m], j,k,k,i, op_collection, lat));
            }

            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_3term, &ta, vec, _1));

            used_elements[m] += 1;
        }

        // 9887 7371 8727

        // 4-fold degenerate (+spin) V_ijil = V_ijli = V_jiil = V_jili  <--- coded
        //                           V_ilij = V_ilji = V_liij = V_liji
        else if ( ((i==k && j!=l) || j==k || (j==l && i!=k)) && (i!=j && k!=l)) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            // add terms for the four possible permutations
            // the terms per permutation correspond to
            // c^dag_{p1, sigma} c^dag_{p2, sigma'} c_{p3, sigma'} d_{p4, sigma}, summed over sigma and sigma'

            term_vec terms;
            append(terms, SSUM::three_term(matrix_elements[m], i,k,l,j, op_collection, lat));
            append(terms, SSUM::three_term(matrix_elements[m], i,l,k,j, op_collection, lat));
            append(terms, SSUM::three_term(matrix_elements[m], j,k,l,i, op_collection, lat));
            append(terms, SSUM::three_term(matrix_elements[m], j,l,k,i, op_collection, lat));

            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_3term, &ta, vec, _1));

            used_elements[m] += 1;
        }

        // 32 (8x4)-fold degenerate V_ijkl = V_jikl = V_ijlk = V_jilk = V_klij = V_lkij = V_klji = V_lkji * spin
        // V_ijkl -> 24 permutations which fall into 3 equivalence classes of 8 elements (with identical V_ matrix element)
        // coded: 4 index permutations which generate all Sz spins
        else if (i!=j && j!=k && k!=l && i!=k && j!=l) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor>  term_vec;

            std::vector<term_descriptor> & vec = this->terms_;

            // add terms for the four possible permutations
            // the terms per permutation correspond to
            // \sum_{sigma, sigma'} c^dag_{p1, sigma} c^dag_{p2, sigma'} c_{p3, sigma'} c_{p4, sigma}

            term_vec terms;
            append(terms, SSUM::four_term(matrix_elements[m], i,k,l,j, op_collection, lat));
            append(terms, SSUM::four_term(matrix_elements[m], i,l,k,j, op_collection, lat));
            append(terms, SSUM::four_term(matrix_elements[m], j,k,l,i, op_collection, lat));
            append(terms, SSUM::four_term(matrix_elements[m], j,l,k,i, op_collection, lat));

            std::for_each(terms.begin(), terms.end(), bind(&ChemHelperSU2<Matrix, SymmGroup>::add_4term, &ta, vec, _1));

            used_elements[m] += 1;
        }

    } // matrix_elements for

    ta.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
