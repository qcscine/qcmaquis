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

#include "dmrg/models/term_descriptor.h"
#include "dmrg/models/chem/util.h"

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

        chem_detail::IndexTuple key = chem_detail::align<S>(i,j,k,l);
        pos_t i_ = key[0], j_ = key[1], k_ = key[2], l_ = key[3];

        if (k_ > l_ && l_ > j_) // eg V_4132
        { // generates up|up|up|up + up|down|down|up + down|up|up|down + down|down|down|down

            ret.push_back(TM::four_term(ops.ident_full.no_couple, 2, -std::sqrt(3.)*matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
            ret.push_back(TM::four_term(ops.ident.no_couple,      1,                matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
        }
        else if (k_ > j_ && j_ > l_) // eg V_4231
        { // generates up|up|up|up + up|down|up|down + down|up|down|up + down|down|down|down

            value_type local_element = matrix_element;
            if (TM::sgn(i,k,l,j)) local_element = -matrix_element;

            ret.push_back(TM::four_term(ops.ident_full.no_couple, 2, std::sqrt(3.)*local_element, i,k,l,j, ops.create, ops.destroy, lat));
            ret.push_back(TM::four_term(ops.ident.no_couple,      1,               local_element, i,k,l,j, ops.create, ops.destroy, lat));
        }
        else if (j_ > k_ && k_ > l_) // eg V_4321
        { // generates up|up|up|up + up|up|down|down + down|down|up|up + down|down|down|down

            ret.push_back(TM::four_term(ops.ident.no_couple, 1, 2.*matrix_element, i,k,l,j, ops.create, ops.destroy, lat));
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
                                   matrix_element, i, j, ops.e2d.no_couple, ops.d2e.no_couple, lat));
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
            false, ops.ident_full.no_couple, std::sqrt(3.) * matrix_element, i, j,
            ops.flip.couple_down, ops.flip.couple_up, ops.flip.couple_down, ops.flip.couple_up, lat
        ));

        ret.push_back(TM::two_term(false, ops.ident.no_couple, -0.5 * matrix_element, i, j, ops.count.no_couple, ops.count.no_couple, lat));

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
                true, ops.ident.no_couple,  std::sqrt(2.)*matrix_element, s, p, ops.create_count.couple_down, ops.create_count.fill_couple_up,
                ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
            ));
        else     // one lonely constructor
            ret.push_back(TM::positional_two_term(
                true, ops.ident.no_couple, -std::sqrt(2.)*matrix_element, s, p, ops.destroy_count.couple_down, ops.destroy_count.fill_couple_up,
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

        ret.push_back(TM::three_term(ops.ident.no_couple, std::sqrt(2.)*matrix_element, same_idx, k, l,
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
            ret.push_back(TM::three_term(ops.ident.no_couple, -std::sqrt(2.)*matrix_element, same_idx, pos1, pos2,
                                         ops.e2d.no_couple, ops.e2d.no_couple,
                                         ops.destroy.couple_down, ops.destroy.fill_couple_up,
                                         ops.destroy.couple_down, ops.destroy.fill_couple_up, lat));
        }

        if (j==l)
        {
            same_idx = j; pos1 = std::min(i, k); pos2 = std::max(i, k);
            ret.push_back(TM::three_term(ops.ident.no_couple, -std::sqrt(2.)*matrix_element, same_idx, pos1, pos2,
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
                    ops.ident_full.no_couple, phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, ops.flip.couple_up, ops.flip.couple_up,
                    ops.create.couple_down, ops.create.fill_couple_down, ops.destroy.couple_down, ops.destroy.fill_couple_down, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, ops.count.no_couple, ops.count.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
            else if (same_idx > std::max(pos1,pos2))
            {
                ret.push_back(TM::three_term(
                    ops.ident_full.no_couple, phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, ops.flip.couple_down, ops.flip.couple_down,
                    ops.create.couple_up, ops.create.fill_couple_up, ops.destroy.couple_up, ops.destroy.fill_couple_up, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, ops.count.no_couple, ops.count.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
            else
            {
                ret.push_back(TM::three_term(
                    ops.ident.no_couple,  phase * std::sqrt(3.)*matrix_element, same_idx, pos1, pos2, ops.flip.no_couple, ops.flip.no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
                ret.push_back(TM::three_term(
                    ops.ident.no_couple, -0.5*std::sqrt(2.)*matrix_element, same_idx, pos1, pos2, ops.count.fill_no_couple, ops.count.fill_no_couple,
                    ops.create.couple_down, ops.create.fill_couple_up, ops.destroy.couple_down, ops.destroy.fill_couple_up, lat
                ));
            }
        }

        return ret;
    }

};

#endif
