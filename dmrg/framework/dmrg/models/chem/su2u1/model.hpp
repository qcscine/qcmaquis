/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

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

    for (subcharge irr=0; irr <= max_irrep; ++irr)
    {
        Index<SymmGroup> phys;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(B, irr), 1));
        phys.insert(std::make_pair(PGCharge<SymmGroup>()(C, irr), 1));
        phys.insert(std::make_pair(D, 1));

        phys_indices.push_back(phys);
    }

    ops = TermMakerSU2<Matrix, SymmGroup>::construct_operators(max_irrep, tag_handler);
    op_collection =  TermMakerSU2<Matrix, SymmGroup>::construct_operator_collection(ops, max_irrep);

}

template <class Matrix, class SymmGroup>
void qc_su2<Matrix, SymmGroup>::create_terms()
{
    typedef typename SymmGroup::subcharge subcharge;
    subcharge N = SymmGroup::particleNumber(this->total_quantum_numbers(parms));

    chem::detail::ChemHelperSU2<Matrix, SymmGroup> ta(parms, lat, tag_handler);
    alps::numeric::matrix<Lattice::pos_t> idx_ = ta.getIdx();
    auto matrix_elements = ta.getMatrixElements();

    std::vector<int> used_elements(matrix_elements.size(), 0);

    typedef TermMakerSU2<Matrix, SymmGroup> TM;

    /**********************************************************************/

    using chem::detail::ChemHelperSU2;
    using chem::detail::append;

    for (std::size_t m=0; m < matrix_elements.size(); ++m) {
        int i = idx_(m, 0);
        int j = idx_(m, 1);
        int k = idx_(m, 2);
        int l = idx_(m, 3);
        auto matrixElement = static_cast<value_type>(matrix_elements[m]);
        // Core electrons energy
        if ( i==-1 && j==-1 && k==-1 && l==-1) {

            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back( std::make_pair(0, ops.ident[lat.get_prop<typename SymmGroup::subcharge>("type", 0)]) );
            this->terms_.push_back(term);

            used_elements[m] += 1;
        }

        // On site energy t_ii
        else if ( i==j && k == -1 && l == -1) {

            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back( std::make_pair(i, ops.count[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
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
                    true, ops.ident, value_type(std::sqrt(2.))*matrixElement,j,i, ops.create, ops.create_fill, ops.destroy, ops.destroy_fill, lat
                ));
                this->terms_.push_back(TermMakerSU2<Matrix, SymmGroup>::positional_two_term(
                    true, ops.ident, value_type(std::sqrt(2.))*matrixElement,i,j, ops.create, ops.create_fill, ops.destroy, ops.destroy_fill, lat
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
                    append(terms, SSUM::three_term(matrixElement * value_type(1./(N-1)), i,kk,kk,j, op_collection, lat));
                    append(terms, SSUM::three_term(matrixElement * value_type(1./(N-1)), j,kk,kk,i, op_collection, lat));
                }

                for (auto&& term: terms)
                    ta.add_3term(vec, term);

                terms.clear();

                append(terms, SSUM::V_term(matrixElement * value_type(1./(N-1)), i,i,i,j, op_collection, lat));
                append(terms, SSUM::V_term(matrixElement * value_type(1./(N-1)), j,i,i,i, op_collection, lat));

                append(terms, SSUM::V_term(matrixElement * value_type(1./(N-1)), i,j,j,j, op_collection, lat));
                append(terms, SSUM::V_term(matrixElement * value_type(1./(N-1)), j,j,j,i, op_collection, lat));

                for (auto&& term: terms)
                    ta.add_2term(vec, term);
            }

            used_elements[m] += 1;
        }

        // On site Coulomb repulsion V_iiii
        else if ( i==j && j==k && k==l) {

            term_descriptor term;
            term.coeff = matrixElement;
            term.push_back(std::make_pair(i, ops.docc[lat.get_prop<typename SymmGroup::subcharge>("type", i)]));
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
            append(terms, SSUM::two_term(matrixElement, s,s,s,p, op_collection, lat));
            append(terms, SSUM::two_term(matrixElement, s,p,s,s, op_collection, lat));

            for (auto&& term: terms)
                ta.add_2term(vec, term);

            used_elements[m] += 1;
        }

        // V_iijj == V_jjii
        else if ( i==j && k==l && j!=k) {

            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            term_vec terms = SSUM::two_term(matrixElement, i,k,k,i, op_collection, lat);

            for (auto&& term: terms)
                ta.add_2term(vec, term);

            used_elements[m] += 1;
        }

        // V_ijij == V_jiji = V_ijji = V_jiij
        else if ( i==k && j==l && i!=j) {
            typedef SpinSumSU2<Matrix, SymmGroup> SSUM;
            typedef std::vector<term_descriptor> term_vec;

            term_vec & vec = this->terms_;

            term_vec terms;
            append(terms, SSUM::two_term(value_type(0.5)*matrixElement, i,i,j,j, op_collection, lat));
            append(terms, SSUM::two_term(value_type(0.5)*matrixElement, j,j,i,i, op_collection, lat));

            append(terms, SSUM::two_term(matrixElement, i,j,i,j, op_collection, lat));

            for (auto&& term: terms)
                ta.add_2term(vec, term);

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
                append(terms, SSUM::three_term(matrixElement, i,k,l,i, op_collection, lat));
                append(terms, SSUM::three_term(matrixElement, i,l,k,i, op_collection, lat));
            }
            else // (k==l)
            {
                append(terms, SSUM::three_term(matrixElement, i,k,k,j, op_collection, lat));
                append(terms, SSUM::three_term(matrixElement, j,k,k,i, op_collection, lat));
            }

            for (auto&& term: terms)
                ta.add_3term(vec, term);

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
            append(terms, SSUM::three_term(matrixElement, i,k,l,j, op_collection, lat));
            append(terms, SSUM::three_term(matrixElement, i,l,k,j, op_collection, lat));
            append(terms, SSUM::three_term(matrixElement, j,k,l,i, op_collection, lat));
            append(terms, SSUM::three_term(matrixElement, j,l,k,i, op_collection, lat));

            for (auto&& term: terms)
                ta.add_3term(vec, term);

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
            append(terms, SSUM::four_term(matrixElement, i,k,l,j, op_collection, lat));
            append(terms, SSUM::four_term(matrixElement, i,l,k,j, op_collection, lat));
            append(terms, SSUM::four_term(matrixElement, j,k,l,i, op_collection, lat));
            append(terms, SSUM::four_term(matrixElement, j,l,k,i, op_collection, lat));

            for (auto&& term: terms) ta.add_4term(vec, term);

            used_elements[m] += 1;
        }

    } // matrix_elements for

    ta.commit_terms(this->terms_);
    maquis::cout << "The hamiltonian will contain " << this->terms_.size() << " terms\n";
}

#endif
