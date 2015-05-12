/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef QC_TERMMAKER_SU2_H
#define QC_TERMMAKER_SU2_H

template <class M, class S>
struct TermMakerSU2 {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef std::vector<tag_type> tag_vec;
    typedef typename term_descriptor::value_type pos_op_t;

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

    static term_descriptor two_term(bool sign, tag_vec full_ident, value_type scale, pos_t i, pos_t j,
                                    tag_vec op1, tag_vec op2, Lattice const & lat)
    {
        term_descriptor term;
        term.is_fermionic = sign;
        term.coeff = scale;
        term.push_back(boost::make_tuple(i, op1[lat.get_prop<sc>("type", i)]));
        term.push_back(boost::make_tuple(j, op2[lat.get_prop<sc>("type", j)]));
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
        term.push_back( boost::make_tuple(start, op1_use[lat.get_prop<sc>("type", start)]) );

        term.push_back( boost::make_tuple(end, op2_use[lat.get_prop<sc>("type", end)]) );

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
        term.push_back( boost::make_tuple(pb, boson_op_use[lat.get_prop<sc>("type", pb)]) );
        term.push_back( boost::make_tuple(start, op1_use[lat.get_prop<sc>("type", start)]) );
        term.push_back( boost::make_tuple(end, op2_use[lat.get_prop<sc>("type", end)]) );

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
            term.push_back(make_tuple(get<0>(sterm[0]), get<1>(sterm[0]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[0]))]));
            term.push_back(make_tuple(get<0>(sterm[1]), get<1>(sterm[1]).couple_up[lat.get_prop<sc>("type", get<0>(sterm[1]))]));
            term.push_back(make_tuple(get<0>(sterm[2]), get<1>(sterm[2]).fill_couple_down[lat.get_prop<sc>("type", get<0>(sterm[2]))]));
            term.push_back(make_tuple(get<0>(sterm[3]), get<1>(sterm[3]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[3]))]));
        }
        else {
            term.push_back(make_tuple(get<0>(sterm[0]), get<1>(sterm[0]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[0]))]));
            term.push_back(make_tuple(get<0>(sterm[1]), get<1>(sterm[1]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[1]))]));
            term.push_back(make_tuple(get<0>(sterm[2]), get<1>(sterm[2]).fill_couple_up[lat.get_prop<sc>("type", get<0>(sterm[2]))]));
            term.push_back(make_tuple(get<0>(sterm[3]), get<1>(sterm[3]).couple_down[lat.get_prop<sc>("type", get<0>(sterm[3]))]));
        }

        return term;
    }
};

template <class M, class S>
struct SpinSumSU2 {

    typedef typename Lattice::pos_t pos_t;
    typedef typename M::value_type value_type;
    typedef ::term_descriptor<value_type> term_descriptor;

    typedef typename TagHandler<M, S>::tag_type tag_type;
    typedef std::vector<tag_type> tag_vec;

    typedef TermMakerSU2<M, S> TM;
    typedef typename TM::OperatorBundle OperatorBundle;
    typedef typename TM::OperatorCollection OperatorCollection;

    //static term_descriptor three_term(tag_vec full_ident,
    //                                  value_type scale, pos_t pb, pos_t p1, pos_t p2,
    //                                  tag_vec boson_op, tag_vec boson_op_fill,
    //                                  tag_vec op1, tag_vec op1_fill, tag_vec op2, tag_vec op2_fill, Lattice const & lat)
    //{
    //    if (i==j || k==l) return three_termA(full_ident, scale, pb, p1, p2, boson_op, boson_op_fill, op1, op1_fill, op2, op2_fill, lat);
    //    else              return three_termB(full_ident, scale, pb, p1, p2, boson_op, boson_op_fill, op1, op1_fill, op2, op2_fill, lat);
    //}

    //static term_descriptor three_termA(tag_vec full_ident,
    //                                   value_type scale, pos_t pb, pos_t p1, pos_t p2,
    //                                   tag_vec boson_op, tag_vec boson_op_fill,
    //                                   tag_vec op1, tag_vec op1_fill, tag_vec op2, tag_vec op2_fill, Lattice const & lat)
    //{
    //}

    static std::vector<term_descriptor> 
    three_term(tag_vec const & ident, tag_vec const & ident_full, value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
               OperatorCollection const & ops, Lattice const & lat)
    {
        tag_vec const & create_couple_up         = ops.create.couple_up        ; 
        tag_vec const & create                   = ops.create.couple_down      ; 
        tag_vec const & create_fill              = ops.create.fill_couple_up   ; 
        tag_vec const & create_fill_couple_down  = ops.create.fill_couple_down ; 

        tag_vec const & destroy_couple_up        = ops.destroy.couple_up       ; 
        tag_vec const & destroy                  = ops.destroy.couple_down     ; 
        tag_vec const & destroy_fill             = ops.destroy.fill_couple_up  ; 
        tag_vec const & destroy_fill_couple_down = ops.destroy.fill_couple_down; 

        tag_vec const & flip_S0                  = ops.flip.no_couple          ; 
        tag_vec const & flip_to_S2               = ops.flip.couple_up          ; 
        tag_vec const & flip_to_S0               = ops.flip.couple_down        ; 

        tag_vec const & count                    = ops.count.no_couple         ; 
        tag_vec const & count_fill               = ops.count.fill_no_couple    ; 

        tag_vec const & e2d                      = ops.e2d.no_couple           ; 
        tag_vec const & d2e                      = ops.d2e.no_couple           ; 

        std::vector<term_descriptor> ret;

        int same_idx, pos1, pos2;

        if (i==k)
        {
            same_idx = i; pos1 = std::min(l, j); pos2 = std::max(l, j);
            ret.push_back(TM::three_term(ident, -std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, e2d, e2d,
            destroy, destroy_fill, destroy, destroy_fill, lat));
        }

        if (j==l)
        {
            same_idx = j; pos1 = std::min(i, k); pos2 = std::max(i, k);
            ret.push_back(TM::three_term(ident, -std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, d2e, d2e,
            create, create_fill, create, create_fill, lat));
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
                    ident_full, phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, flip_to_S2, flip_to_S2,
                    create, create_fill_couple_down, destroy, destroy_fill_couple_down, lat
                ));
                ret.push_back(TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, count, count,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
            else if (same_idx > std::max(pos1,pos2))
            {
                ret.push_back(TM::three_term(
                    ident_full, phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, flip_to_S0, flip_to_S0,
                    create_couple_up, create_fill, destroy_couple_up, destroy_fill, lat
                ));
                ret.push_back(TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, count, count, create, create_fill, destroy, destroy_fill, lat
                ));
            }
            else
            {
                ret.push_back(TM::three_term(
                    ident,  phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, flip_S0, flip_S0, create, create_fill, destroy, destroy_fill, lat
                ));
                ret.push_back(TM::three_term(
                    ident, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, count_fill, count_fill, create, create_fill, destroy, destroy_fill, lat
                ));
            }
        }

        return ret;
    }

    static std::vector<term_descriptor> 
    four_term(tag_vec const & ident, tag_vec const & ident_full, value_type matrix_element, pos_t i, pos_t k, pos_t l, pos_t j,
              OperatorBundle const & create_pkg, OperatorBundle const & destroy_pkg, Lattice const & lat)
    {
        std::vector<term_descriptor> ret;

        // These 3 cases produce different S_z spin patterns, which differ along with different index permutations
        // As in standard notation of the Hamiltonian, the first two positions get a creator, the last two a destructor

        chem_detail::IndexTuple key = chem_detail::align<S>(i,j,k,l);
        pos_t i_ = key[0], j_ = key[1], k_ = key[2], l_ = key[3];

        if (k_ > l_ && l_ > j_) // eg V_4132
        { // generates up|up|up|up + up|down|down|up + down|up|up|down + down|down|down|down

            ret.push_back(TM::four_term(ident_full, 2, -std::sqrt(3.)*matrix_element, i,k,l,j, create_pkg, destroy_pkg, lat));
            ret.push_back(TM::four_term(ident,      1,                matrix_element, i,k,l,j, create_pkg, destroy_pkg, lat));
        }
        else if (k_ > j_ && j_ > l_) // eg V_4231
        { // generates up|up|up|up + up|down|up|down + down|up|down|up + down|down|down|down

            value_type local_element = matrix_element;
            if (TM::sgn(i,k,l,j)) local_element = -matrix_element;

            ret.push_back(TM::four_term(ident_full, 2, std::sqrt(3.)*local_element, i,k,l,j, create_pkg, destroy_pkg, lat));
            ret.push_back(TM::four_term(ident,      1,               local_element, i,k,l,j, create_pkg, destroy_pkg, lat));
        }
        else if (j_ > k_ && k_ > l_) // eg V_4321
        { // generates up|up|up|up + up|up|down|down + down|down|up|up + down|down|down|down

            ret.push_back(TM::four_term(ident, 1, 2.*matrix_element, i,k,l,j, create_pkg, destroy_pkg, lat));
        }
        else { throw std::runtime_error("unexpected index arrangment in V_ijkl term\n"); }

        return ret;
    }

};

#endif
