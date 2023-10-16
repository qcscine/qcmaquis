/**
 * @file
 * @copyright This code is licensed under the 3-clause BSD license.
 *            Copyright ETH Zurich, Laboratory of Physical Chemistry, Reiher Group.
 *            See LICENSE.txt for details.
 */

#ifndef QC_TERMMAKER_SU2_H
#define QC_TERMMAKER_SU2_H

#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/chem/util.h"

template <class M, class S>
struct TermMakerSU2 {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;
    typedef typename operator_selector<M, S>::type op_t;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef std::vector<tag_type> tag_vec;
    typedef typename term_descriptor::value_type pos_op_t;
    typedef std::shared_ptr<TagHandler<M, S> > tag_handler_t;

    typedef typename S::subcharge sc;

    struct OperatorBundle
    {
        tag_vec couple_up;
        tag_vec couple_down;
        tag_vec fill_couple_up;
        tag_vec fill_couple_down;

        tag_vec no_couple;
        tag_vec fill_no_couple;
    };

    // Operators replaces loose tag_vecs in qc_su2
    struct Operators
    {
        tag_vec create_fill,
                create,
                destroy_fill,
                destroy,
                create_fill_couple_down,
                destroy_fill_couple_down,
                create_couple_up,
                destroy_couple_up,
                create_fill_count,
                create_count,
                destroy_fill_count,
                destroy_count,
                count,
                docc,
                e2d,
                d2e,
                flip_S0,
                flip_to_S2,
                flip_to_S0,
                ident,
                ident_full,
                fill,
                count_fill;
    };

    struct OperatorCollection
    {
        OperatorBundle ident;
        OperatorBundle ident_full;
        OperatorBundle fill;
        OperatorBundle create;
        OperatorBundle destroy;
        OperatorBundle count;
        OperatorBundle flip;
        OperatorBundle e2d;
        OperatorBundle d2e;
        OperatorBundle docc;
        OperatorBundle create_count;
        OperatorBundle destroy_count;
    };

    typedef boost::tuple<pos_t, OperatorBundle> pos_bundle_t;

    template <class Tuple>
    static bool compare_tag(Tuple p1,
                            Tuple p2)
    {
        return boost::tuples::get<0>(p1) < boost::tuples::get<0>(p2);
    }

    template <class I>
    static bool sgn(I i, I j, I k, I l)
    {
        // Simple O(n^2) algorithm to determine sign of permutation
        I idx[] = { i,j,k,l };
        I inv_count=0, n=4;
        for(I c1 = 0; c1 < n - 1; c1++)
            for(I c2 = c1+1; c2 < n; c2++)
                if(idx[c1] > idx[c2]) inv_count++;

        return (inv_count % 2 != 0);
    }

    // Functions needed for generating operators in qc_su2
    static std::vector<op_t> generate_site_specific_ops(op_t const & op, typename S::subcharge max_irrep)
    {
        PGDecorator<S> set_symm;
        std::vector<op_t> ret;
        for (typename S::subcharge sc=0; sc < max_irrep+1; ++sc) {
            op_t mod(set_symm(op.basis(), sc));
            mod.spin() = op.spin();
            for (std::size_t b = 0; b < op.n_blocks(); ++b)
                mod[b] = op[b];

            ret.push_back(mod);
        }
        return ret;
    }

    static tag_vec register_site_specific(std::vector<op_t> const & ops, tag_detail::operator_kind kind, typename S::subcharge max_irrep, const tag_handler_t & tag_handler)
    {
        std::vector<tag_type> ret;
        for (typename S::subcharge sc=0; sc < max_irrep+1; ++sc) {
            std::pair<tag_type, value_type> newtag = tag_handler->checked_register(ops[sc], kind);
            assert( newtag.first < tag_handler->size() );
            assert( std::abs(newtag.second - value_type(1.)) == value_type() );
            ret.push_back(newtag.first);
        }

        return ret;
    }

    static Operators construct_operators(typename S::subcharge max_irrep, const tag_handler_t & tag_handler)
    {
        Operators ret;

        typename S::charge A(0), B(0), C(0), D(0);
        A[0] = 2; // 20
        B[0] = 1; B[1] =  1; // 11
        C[0] = 1; C[1] = -1; // 1-1
        // D = 00

        SpinDescriptor<symm_traits::SU2Tag> one_half_up(1,0,1);
        SpinDescriptor<symm_traits::SU2Tag> one_half_down(1,1,0);
        SpinDescriptor<symm_traits::SU2Tag> one_up(2,0,2);
        SpinDescriptor<symm_traits::SU2Tag> one_flat(2,1,1);
        SpinDescriptor<symm_traits::SU2Tag> one_down(2,2,0);

        // cheaper to use this for spin0 tensors, instead of ident_full
        op_t ident_op;
        ident_op.insert_block(M(1,1,1), A, A);
        ident_op.insert_block(M(1,1,1), B, B);
        ident_op.insert_block(M(1,1,1), C, C);
        ident_op.insert_block(M(1,1,1), D, D);

        // apply if spin > 0
        op_t ident_full_op;
        ident_full_op.insert_block(M(1,1,1), A, A);
        ident_full_op.insert_block(M(1,1,1), D, D);
        ident_full_op.insert_block(M(1,1,1), B, B);
        ident_full_op.insert_block(M(1,1,1), C, C);
        ident_full_op.insert_block(M(1,1,1), B, C);
        ident_full_op.insert_block(M(1,1,1), C, B);

        op_t fill_op;
        fill_op.insert_block(M(1,1,1),  A, A);
        fill_op.insert_block(M(1,1,1),  D, D);
        fill_op.insert_block(M(1,1,-1), B, B);
        fill_op.insert_block(M(1,1,-1), C, C);
        fill_op.insert_block(M(1,1,-1), B, C);
        fill_op.insert_block(M(1,1,-1), C, B);

        /*************************************************************/

        op_t create_fill_op;
        create_fill_op.spin() = one_half_up;
        create_fill_op.insert_block(M(1,1,sqrt(2.)), B, A);
        create_fill_op.insert_block(M(1,1,sqrt(2.)), C, A);
        create_fill_op.insert_block(M(1,1,1), D, B);
        create_fill_op.insert_block(M(1,1,1), D, C);

        op_t destroy_op;
        destroy_op.spin() = one_half_down;
        destroy_op.insert_block(M(1,1,1), A, B);
        destroy_op.insert_block(M(1,1,1), A, C);
        destroy_op.insert_block(M(1,1,sqrt(2.)), B, D);
        destroy_op.insert_block(M(1,1,sqrt(2.)), C, D);

        op_t destroy_fill_op;
        destroy_fill_op.spin() = one_half_up;
        destroy_fill_op.insert_block(M(1,1,1), A, B);
        destroy_fill_op.insert_block(M(1,1,1), A, C);
        destroy_fill_op.insert_block(M(1,1,-sqrt(2.)), B, D);
        destroy_fill_op.insert_block(M(1,1,-sqrt(2.)), C, D);

        op_t create_op;
        create_op.spin() = one_half_down;
        create_op.insert_block(M(1,1,sqrt(2.)), B, A);
        create_op.insert_block(M(1,1,sqrt(2.)), C, A);
        create_op.insert_block(M(1,1,-1), D, B);
        create_op.insert_block(M(1,1,-1), D, C);

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
        create_fill_count_op.insert_block(M(1,1,sqrt(2.)), B, A);
        create_fill_count_op.insert_block(M(1,1,sqrt(2.)), C, A);

        op_t destroy_count_op;
        destroy_count_op.spin() = one_half_down;
        destroy_count_op.insert_block(M(1,1,1), A, B);
        destroy_count_op.insert_block(M(1,1,1), A, C);

        op_t destroy_fill_count_op;
        destroy_fill_count_op.spin() = one_half_up;
        destroy_fill_count_op.insert_block(M(1,1,1), A, B);
        destroy_fill_count_op.insert_block(M(1,1,1), A, C);

        op_t create_count_op;
        create_count_op.spin() = one_half_down;
        create_count_op.insert_block(M(1,1,sqrt(2.)), B, A);
        create_count_op.insert_block(M(1,1,sqrt(2.)), C, A);

        /*************************************************************/

        op_t count_op;
        count_op.insert_block(M(1,1,2), A, A);
        count_op.insert_block(M(1,1,1), B, B);
        count_op.insert_block(M(1,1,1), C, C);

        op_t docc_op;
        docc_op.insert_block(M(1,1,1), A, A);

        op_t e2d_op;
        e2d_op.insert_block(M(1,1,1), D, A);

        op_t d2e_op;
        d2e_op.insert_block(M(1,1,1), A, D);

        op_t count_fill_op;
        count_fill_op.insert_block(M(1,1,2),  A, A);
        count_fill_op.insert_block(M(1,1,-1), B, B);
        count_fill_op.insert_block(M(1,1,-1), C, C);
        count_fill_op.insert_block(M(1,1,-1), B, C);
        count_fill_op.insert_block(M(1,1,-1), C, B);

        op_t flip_to_S2_op;
        flip_to_S2_op.spin() = one_up;
        flip_to_S2_op.insert_block(M(1,1,std::sqrt(3./2)), B, B);
        flip_to_S2_op.insert_block(M(1,1,std::sqrt(3./2.)), C, C);
        flip_to_S2_op.insert_block(M(1,1,std::sqrt(3./2.)),  B, C);
        flip_to_S2_op.insert_block(M(1,1,std::sqrt(3./2.)),  C, B);

        op_t flip_to_S0_op = flip_to_S2_op;
        flip_to_S0_op.spin() = one_down;

        op_t flip_S0_op = flip_to_S2_op;
        flip_S0_op.spin() = one_flat;

        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/

        #define GENERATE_SITE_SPECIFIC(opname) std::vector<op_t> opname ## s = generate_site_specific_ops(opname, max_irrep);

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

        #undef GENERATE_SITE_SPECIFIC

        #define REGISTER(op, kind) ret.op = register_site_specific(op ## _ops, kind, max_irrep, tag_handler);

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

        #define HERMITIAN(op1, op2) for (int hh=0; hh < ret.op1.size(); ++hh) tag_handler->hermitian_pair(ret.op1[hh], ret.op2[hh]);
        HERMITIAN(create_fill, destroy_fill)
        HERMITIAN(create, destroy)
        HERMITIAN(e2d, d2e)

        HERMITIAN(create_fill_count, destroy_fill_count) // useless
        HERMITIAN(create_count, destroy_count)

        HERMITIAN(create_fill_couple_down, destroy_fill_couple_down) // useless

        HERMITIAN(create_couple_up, destroy_couple_up)
        #undef HERMITIAN

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

        return ret;
    }

    static OperatorCollection construct_operator_collection(Operators const& ops, typename S::subcharge max_irrep)
    {
        OperatorCollection op_collection;

        OperatorBundle create_pkg, destroy_pkg;
        OperatorBundle create_count_pkg, destroy_count_pkg;

        create_pkg.couple_up = ops.create_couple_up;
        create_pkg.couple_down = ops.create;
        create_pkg.fill_couple_up = ops.create_fill;
        create_pkg.fill_couple_down = ops.create_fill_couple_down;

        destroy_pkg.couple_up = ops.destroy_couple_up;
        destroy_pkg.couple_down = ops.destroy;
        destroy_pkg.fill_couple_up = ops.destroy_fill;
        destroy_pkg.fill_couple_down = ops.destroy_fill_couple_down;

        create_count_pkg.couple_down = ops.create_count;
        create_count_pkg.fill_couple_up = ops.create_fill_count;

        destroy_count_pkg.couple_down = ops.destroy_count;
        destroy_count_pkg.fill_couple_up = ops.destroy_fill_count;

        /**********************************************************************/

        op_collection.ident     .no_couple = ops.ident;
        op_collection.ident_full.no_couple = ops.ident_full;
        op_collection.fill      .no_couple = ops.fill;

        op_collection.create               = create_pkg;
        op_collection.destroy              = destroy_pkg;

        op_collection.count     .no_couple = ops.count;
        op_collection.count     .fill_no_couple = ops.count_fill;

        op_collection.create_count         = create_count_pkg;
        op_collection.destroy_count        = destroy_count_pkg;

        op_collection.e2d       .no_couple = ops.e2d;
        op_collection.d2e       .no_couple = ops.d2e;
        op_collection.docc      .no_couple = ops.docc;

        op_collection.flip      .no_couple = ops.flip_S0;
        op_collection.flip      .couple_up = ops.flip_to_S2;
        op_collection.flip      .couple_down = ops.flip_to_S0;

        /**********************************************************************/
        return op_collection;
    }

    static term_descriptor two_term(bool sign, tag_vec full_ident, value_type scale, pos_t i, pos_t j,
                                    tag_vec op1, tag_vec op2, Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;
        term.push_back(std::make_pair(i, op1[lat.get_prop<sc>("type", i)]));
        term.push_back(std::make_pair(j, op2[lat.get_prop<sc>("type", j)]));
        return term;
    }

    static term_descriptor positional_two_term(bool sign, tag_vec full_ident, value_type scale, pos_t i, pos_t j,
                                     tag_vec op1, tag_vec op1_fill, tag_vec op2, tag_vec op2_fill, Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;

        tag_vec op1_use = (i<j) ? op1_fill : op2_fill;
        tag_vec op2_use = (i<j) ? op2 : op1;
        if (j<i && sign) term.coeff = -term.coeff;

        pos_t start = std::min(i,j), end = std::max(i,j);
        term.push_back( std::make_pair(start, op1_use[lat.get_prop<sc>("type", start)]) );

        term.push_back( std::make_pair(end, op2_use[lat.get_prop<sc>("type", end)]) );

        return term;
    }

    static term_descriptor three_term(tag_vec full_ident,
                                      value_type scale, pos_t pb, pos_t p1, pos_t p2,
                                      tag_vec boson_op, tag_vec boson_op_fill,
                                      tag_vec op1, tag_vec op1_fill, tag_vec op2, tag_vec op2_fill, Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = true;
        term.coeff = scale;

        tag_vec boson_op_use, op1_use, op2_use;

        op1_use = (p1<p2) ? op1_fill : op2_fill;
        op2_use = (p1<p2) ? op2 : op1;

        if ( (pb>p1 && pb<p2) || (pb>p2 && pb<p1) ) {
            // if the bosonic operator is in between
            // the fermionic operators, use fill-multiplied version
            boson_op_use = boson_op_fill;
        }
        else {
            boson_op_use = boson_op;
        }

        if (p2<p1) term.coeff = -term.coeff;

        pos_t start = std::min(p1,p2), end = std::max(p1,p2);
        term.push_back( std::make_pair(pb, boson_op_use[lat.get_prop<sc>("type", pb)]) );
        term.push_back( std::make_pair(start, op1_use[lat.get_prop<sc>("type", start)]) );
        term.push_back( std::make_pair(end, op2_use[lat.get_prop<sc>("type", end)]) );

        //maquis::cout << "  term " << pb << p1 << p2 << "   " << boson_op_use[lat.get_prop<sc>("type", pb)] << "," << op1_use[lat.get_prop<sc>("type", start)] << "," << op2_use[lat.get_prop<sc>("type", end)] << std::endl;

        // sort the terms w.r. to position
        term.canonical_order();

        return term;
    }

    static term_descriptor four_term(tag_vec full_ident, int max_two_S,
                                     value_type scale, pos_t i, pos_t j, pos_t k, pos_t l,
                                     OperatorBundle op_i, OperatorBundle op_k, Lattice const & lat)
    {
        using boost::tuples::get;
        using boost::make_tuple;

        term_descriptor term;
        term.is_fermionic = true;
        term.coeff = scale;

        if (sgn(i,j,k,l))
            term.coeff = -term.coeff;

        std::vector<pos_bundle_t> sterm;
        sterm.push_back(make_tuple(i, op_i));
        sterm.push_back(make_tuple(j, op_i));
        sterm.push_back(make_tuple(k, op_k));
        sterm.push_back(make_tuple(l, op_k));
        std::sort(sterm.begin(), sterm.end(), compare_tag<pos_bundle_t>);

        if (max_two_S == 2) {
            term.push_back(std::make_pair(get<0>(sterm[0]), get<1>(sterm[0]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[0]))]));
            term.push_back(std::make_pair(get<0>(sterm[1]), get<1>(sterm[1]).couple_up[lat.get_prop<sc>("type", get<0>(sterm[1]))]));
            term.push_back(std::make_pair(get<0>(sterm[2]), get<1>(sterm[2]).fill_couple_down[lat.get_prop<sc>("type", get<0>(sterm[2]))]));
            term.push_back(std::make_pair(get<0>(sterm[3]), get<1>(sterm[3]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[3]))]));
        }
        else {
            term.push_back(std::make_pair(get<0>(sterm[0]), get<1>(sterm[0]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[0]))]));
            term.push_back(std::make_pair(get<0>(sterm[1]), get<1>(sterm[1]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[1]))]));
            term.push_back(std::make_pair(get<0>(sterm[2]), get<1>(sterm[2]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[2]))]));
            term.push_back(std::make_pair(get<0>(sterm[3]), get<1>(sterm[3]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[3]))]));
        }

        return term;
    }
};

template <class M, class S>
class SpinSumSU2 {

public:
    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef std::vector<tag_type> tag_vec;

    typedef TermMakerSU2<M, S> TM;
    typedef typename TM::OperatorBundle OperatorBundle;
    typedef typename TM::OperatorCollection OperatorCollection;

    static std::vector<term_descriptor>
    V_term(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
         OperatorCollection const & ops, Lattice const & lat)
    {
        Lattice::pos_t p[4] = {i,k,l,j};
        std::vector<pos_t> ps(p, p + 4);
        std::sort(ps.begin(), ps.end());
        std::vector<pos_t>::iterator it = std::unique(ps.begin(), ps.end());

        std::size_t n_unique = std::distance(ps.begin(), it);

        switch(n_unique) {
            case 4:
                return four_term(matrix_element, i, k, l, j, ops, lat);
            case 3:
                return three_term(matrix_element, i, k, l, j, ops, lat);
            case 2:
                return two_term(matrix_element, i, k, l, j, ops, lat);
            case 1:
                term_descriptor term;
                term.coeff = value_type(2.)*matrix_element; // 2 spin combinations are non-zero
                term.push_back(std::make_pair(i, ops.docc.no_couple[lat.get_prop<typename S::subcharge>("type", i)]));
                return std::vector<term_descriptor>(1, term);
        }

        return std::vector<term_descriptor>();
    }

    static std::vector<term_descriptor>
    two_term(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
             OperatorCollection const & ops, Lattice const & lat)
    {
        if (i==j && k==l && j!=k) return two_termA(matrix_element, i, k, l, j, ops, lat);
        else if (i==k && j==l && j!=k) return two_termB1(matrix_element, i, k, l, j, ops, lat);
        else if (i==l && j==k && i!=j) return two_termB2(matrix_element, i, k, l, j, ops, lat);
        else if ((i==k && k==l) ||
                 (k==l && l==j) ||
                 (i==l && l==j) ||
                 (i==k && k==j)) return two_termC(matrix_element, i, k, l, j, ops, lat);
        else { throw std::runtime_error("Unexpected index arrangement for V_ijjj term\n"); }
    }

    static std::vector<term_descriptor>
    three_term(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
               OperatorCollection const & ops, Lattice const & lat)
    {
        if (i==j || k==l) return three_termA(matrix_element, i, k, l, j, ops, lat);
        else return three_termB(matrix_element, i, k, l, j, ops, lat);
    }

    static std::vector<term_descriptor>
    four_term(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;

        // These 3 cases produce different S_z spin patterns, which differ along with different index permutations
        // As in standard notation of the Hamiltonian, the first two positions get a creator, the last two a destructor

        chem::detail::IndexTuple key = chem::detail::align<S>(i,j,k,l);
        pos_t i_ = key[0], j_ = key[1], k_ = key[2], l_ = key[3];

        if (k_ > l_ && l_ > j_) // eg V_4132
        { // generates up|up|up|up + up|down|down|up + down|up|up|down + down|down|down|down

            ret.push_back(TM::four_term(ops.ident_full.no_couple, 2, value_type(-std::sqrt(3.))*matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
            ret.push_back(TM::four_term(ops.ident.no_couple,      1,                matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
        }
        else if (k_ > j_ && j_ > l_) // eg V_4231
        { // generates up|up|up|up + up|down|up|down + down|up|down|up + down|down|down|down

            value_type local_element = matrix_element;
            if (TM::sgn(i,k,l,j)) local_element = -matrix_element;

            ret.push_back(TM::four_term(ops.ident_full.no_couple, 2, value_type(std::sqrt(3.))*local_element, i,k,l,j, ops.create, ops.destroy, lat));
            ret.push_back(TM::four_term(ops.ident.no_couple,      1,               local_element, i,k,l,j, ops.create, ops.destroy, lat));
        }
        else if (j_ > k_ && k_ > l_) // eg V_4321
        { // generates up|up|up|up + up|up|down|down + down|down|up|up + down|down|down|down

            ret.push_back(TM::four_term(ops.ident.no_couple, 1, value_type(2.)*matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
        }
        else { throw std::runtime_error("unexpected index arrangment in V_ijkl term\n"); }

        return ret;
    }

private:

    static std::vector<term_descriptor>
    two_termA(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;
        ret.push_back(TM::two_term(false, ops.ident.no_couple,
                                   matrix_element, i, k, ops.count.no_couple, ops.count.no_couple, lat));
        return ret;
    }

    static std::vector<term_descriptor>
    two_termB1(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;
        ret.push_back(TM::two_term(false, ops.ident.no_couple,
                                   value_type(2.0) * matrix_element, i, j, ops.e2d.no_couple, ops.d2e.no_couple, lat));
        return ret;
    }

    static std::vector<term_descriptor>
    two_termB2(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;

        // here we have spin0--j--spin1--i--spin0
        // the sqrt(3.) counteracts the Clebsch coeff C^{110}_{mrm'} which applies when the spin1 couples back to spin0
        ret.push_back(TM::positional_two_term(
            false, ops.ident_full.no_couple, value_type(std::sqrt(3.)) * matrix_element, i, j,
            ops.flip.couple_down, ops.flip.couple_up, ops.flip.couple_down, ops.flip.couple_up, lat
        ));

        ret.push_back(TM::two_term(false, ops.ident.no_couple, value_type(-0.5) * matrix_element, i, j, ops.count.no_couple, ops.count.no_couple, lat));

        return ret;
    }

    static std::vector<term_descriptor>
    two_termC(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorCollection const & ops, Lattice const & lat)
    {
        int s, p;

        if (i==k && k==l)      { s = i; p = j; }
        else if (k==l && l==j) { s = j; p = i; }
        else if (i==l && l==j) { s = i; p = k; }
        else if (i==k && k==j) { s = i; p = l; }
        else  { throw std::runtime_error("Term generation logic has failed for V_ijjj term\n"); }

        std::vector<term_descriptor> ret;

        if (i==k) // one lonely destructor
            ret.push_back(TM::positional_two_term(
                true, ops.ident.no_couple,  value_type(std::sqrt(2.))*matrix_element, s, p, ops.create_count.couple_down, ops.create_count.fill_couple_up,
                ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
            ));
        else     // one lonely constructor
            ret.push_back(TM::positional_two_term(
                true, ops.ident.no_couple, value_type(-std::sqrt(2.))*matrix_element, s, p, ops.destroy_count.couple_down, ops.destroy_count.fill_couple_up,
                ops.create.couple_down, ops.create.fill_couple_up, lat
            ));

        return ret;
    }

    static std::vector<term_descriptor>
    three_termA(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
               OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;
        int same_idx;
        if (i==j) { same_idx = i; }
        if (k==l) { same_idx = k; k = i; l = j; }

        ret.push_back(TM::three_term(ops.ident.no_couple, value_type(std::sqrt(2.))*matrix_element, same_idx, k, l,
                                     ops.count.no_couple, ops.count.fill_no_couple,
                                     ops.create.couple_down, ops.create.fill_couple_up,
                                     ops.destroy.couple_down, ops.destroy.fill_couple_up, lat));
        return ret;
    }

    static std::vector<term_descriptor>
    three_termB(value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
                OperatorCollection const & ops, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;

        int same_idx, pos1, pos2;

        if (i==k)
        {
            same_idx = i; pos1 = std::min(l, j); pos2 = std::max(l, j);
            ret.push_back(TM::three_term(ops.ident.no_couple, value_type(-std::sqrt(2.))*matrix_element, same_idx, pos1, pos2,
                                         ops.e2d.no_couple, ops.e2d.no_couple,
                                         ops.destroy.couple_down, ops.destroy.fill_couple_up,
                                         ops.destroy.couple_down, ops.destroy.fill_couple_up, lat));
        }

        if (j==l)
        {
            same_idx = j; pos1 = std::min(i, k); pos2 = std::max(i, k);
            ret.push_back(TM::three_term(ops.ident.no_couple, value_type(-std::sqrt(2.))*matrix_element, same_idx, pos1, pos2,
                                         ops.d2e.no_couple, ops.d2e.no_couple,
                                         ops.create.couple_down, ops.create.fill_couple_up,
                                         ops.create.couple_down, ops.create.fill_couple_up, lat));
        }

        if (j==k || i==l)
        {
            if (j==k) { same_idx = j; pos1 = l; pos2 = i; }
            else { same_idx = i; pos1 = j; pos2 = k; }

            value_type phase = 1.;
            if(TM::sgn(i,k,l,j))
                phase = -1.;

            if ( same_idx < std::min(pos1,pos2) )
            {
                ret.push_back(TM::three_term(
                    ops.ident_full.no_couple, phase * value_type(std::sqrt(3.))*matrix_element, same_idx, pos1, pos2, ops.flip.couple_up, ops.flip.couple_up,
                    ops.create.couple_down, ops.create.fill_couple_down, ops.destroy.couple_down, ops.destroy.fill_couple_down, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, value_type(-0.5*std::sqrt(2.))*matrix_element, same_idx, pos1, pos2, ops.count.no_couple, ops.count.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
            else if (same_idx > std::max(pos1,pos2))
            {
                ret.push_back(TM::three_term(
                    ops.ident_full.no_couple, phase * value_type(std::sqrt(3.))*matrix_element, same_idx, pos1, pos2, ops.flip.couple_down, ops.flip.couple_down,
                    ops.create.couple_up, ops.create.fill_couple_up, ops.destroy.couple_up, ops.destroy.fill_couple_up, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, value_type(-0.5*std::sqrt(2.))*matrix_element, same_idx, pos1, pos2, ops.count.no_couple, ops.count.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
            else
            {
                ret.push_back(TM::three_term(
                    ops.ident.no_couple,  phase * value_type(std::sqrt(3.))*matrix_element, same_idx, pos1, pos2, ops.flip.no_couple, ops.flip.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, value_type(-0.5*std::sqrt(2.))*matrix_element, same_idx, pos1, pos2, ops.count.fill_no_couple, ops.count.fill_no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
        }

        return ret;
    }

};

#endif
