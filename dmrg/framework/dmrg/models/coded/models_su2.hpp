/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
 *               2014-2014 by Sebastian Keller <sebkelle@phys.ethz.ch>
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

#ifndef MODELS_CODED_SU2_H
#define MODELS_CODED_SU2_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** FERMI HUBBARD */
template<class Matrix>
class FermiHubbardSU2 : public model_impl<Matrix, SU2U1>
{
public:
    typedef model_impl<Matrix, SU2U1> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;

    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Matrix::value_type value_type;
    typedef Lattice::pos_t pos_t;
    
    FermiHubbardSU2(const Lattice& lat_, BaseParameters & parms_)
    : lat(lat_)
    , parms(parms_)
    , tag_handler(new TagHandler<Matrix, SU2U1>())
    {
        SU2U1::charge A(0), B(0), C(0), D(0);
        A[0] = 2; // 20
        B[0] = 1; B[1] =  1; // 11
        C[0] = 1; C[1] = -1; // 1-1
        // D = 00

        SpinDescriptor<typename symm_traits::SymmType<SU2U1>::type> one_half_up(1,1);
        SpinDescriptor<typename symm_traits::SymmType<SU2U1>::type> one_half_down(1,-1);

        phys.insert(std::make_pair(A,1));
        phys.insert(std::make_pair(B,1));
        phys.insert(std::make_pair(C,1));
        phys.insert(std::make_pair(D,1));

        op_t identity_op;
        identity_op.insert_block(Matrix(1,1,1), A, A);
        identity_op.insert_block(Matrix(1,1,1), B, B);
        identity_op.insert_block(Matrix(1,1,1), C, C);
        identity_op.insert_block(Matrix(1,1,1), D, D);

        op_t fill_op;
        fill_op.insert_block(Matrix(1,1,1),  A, A);
        fill_op.insert_block(Matrix(1,1,1),  D, D);
        fill_op.insert_block(Matrix(1,1,-1), B, B); 
        fill_op.insert_block(Matrix(1,1,-1), C, C); 
        fill_op.insert_block(Matrix(1,1,-1),  B, C); 
        fill_op.insert_block(Matrix(1,1,-1),  C, B); 

        /*************************************************************/

        op_t create_fill_op;
        create_fill_op.spin = one_half_up;
        create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_fill_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
        create_fill_op.insert_block(Matrix(1,1,1), D, B);
        create_fill_op.insert_block(Matrix(1,1,1), D, C);

        op_t destroy_op;
        destroy_op.spin = one_half_down;
        destroy_op.insert_block(Matrix(1,1,1), A, B);
        destroy_op.insert_block(Matrix(1,1,1), A, C);
        destroy_op.insert_block(Matrix(1,1,sqrt(2.)), B, D);
        destroy_op.insert_block(Matrix(1,1,sqrt(2.)), C, D);

        op_t destroy_fill_op;
        destroy_fill_op.spin = one_half_up;
        destroy_fill_op.insert_block(Matrix(1,1,1), A, B);
        destroy_fill_op.insert_block(Matrix(1,1,1), A, C);
        destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), B, D);
        destroy_fill_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

        op_t create_op;
        create_op.spin = one_half_down;
        create_op.insert_block(Matrix(1,1,sqrt(2.)), B, A);
        create_op.insert_block(Matrix(1,1,sqrt(2.)), C, A);
        create_op.insert_block(Matrix(1,1,-1), D, B);
        create_op.insert_block(Matrix(1,1,-1), D, C);

        /*************************************************************/

        op_t count_op;
        count_op.insert_block(Matrix(1,1,2), A, A);
        count_op.insert_block(Matrix(1,1,1), B, B);
        count_op.insert_block(Matrix(1,1,1), C, C);

        op_t doubly_occ_op;
        doubly_occ_op.insert_block(Matrix(1,1,1), A, A);

        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(identity,     tag_detail::bosonic)
        REGISTER(fill,         tag_detail::bosonic)
        REGISTER(create_fill,  tag_detail::fermionic)
        REGISTER(create,       tag_detail::fermionic)
        REGISTER(destroy,      tag_detail::fermionic)
        REGISTER(destroy_fill, tag_detail::fermionic)
        REGISTER(count,        tag_detail::bosonic)
        REGISTER(doubly_occ,   tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        struct TM {
            static term_descriptor positional_two_term(bool sign, tag_type full_ident, value_type scale, pos_t i, pos_t j,
                                                       tag_type op1, tag_type op1_fill, tag_type op2, tag_type op2_fill)
            {
                term_descriptor term;
                term.is_fermionic = sign;
                term.full_identity = full_ident;
                term.coeff = scale;

                tag_type op1_use = (i<j) ? op1_fill : op2_fill;
                tag_type op2_use = (i<j) ? op2 : op1;
                if (j<i && sign) term.coeff = -term.coeff;

                int start = std::min(i,j), end = std::max(i,j);
                term.push_back( boost::make_tuple(start, op1_use) );

                term.push_back( boost::make_tuple(end, op2_use) );

                return term;
            }
        };
        /**********************************************************************/
        
        value_type U = parms["U"];
        std::pair<tag_type, value_type> ptag;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                term_descriptor term;
                term.coeff = U;
                term.push_back( boost::make_tuple(p, doubly_occ) );
                this->terms_.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                value_type ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));

                // The sqrt(2.) balances the magnitudes of Clebsch coeffs C^{1/2 1/2 0}_{mrm'} which apply at the second spin-1/2 operator
                this->terms_.push_back(TM::positional_two_term(
                    true, identity, std::sqrt(2.) * ti, *hopto, p, create, create_fill, destroy, destroy_fill
                ));
                this->terms_.push_back(TM::positional_two_term(
                    true, identity, std::sqrt(2.) * ti, p, *hopto, create, create_fill, destroy, destroy_fill
                ));
            }
        }
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<SU2U1> const & phys_dim(size_t type) const
    {
        return phys;
    }
    
    measurements_type measurements () const
    {
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;

        /*
        if (parms["ENABLE_MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, SU2U1>("Density up", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_up))) );
        }
        if (parms["ENABLE_MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, SU2U1>("Density down", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_down))) );
        }
        
        if (parms["ENABLE_MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, SU2U1>("Local density up", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_up))) );
        }
        if (parms["ENABLE_MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, SU2U1>("Local density down", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_down))) );
        }
        
        if (parms["ENABLE_MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create_up)), true) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy_up)), true) );
            meas.push_back( new measurements::correlations<Matrix, SU2U1>("Onebody density matrix up", lat,
                                                                       op_vec(1,this->identity_matrix(0)),
                                                                       op_vec(1,this->filling_matrix(0)),
                                                                       ops, false, false) );
        }
        if (parms["ENABLE_MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create_down)), true) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy_down)), true) );
            meas.push_back( new measurements::correlations<Matrix, SU2U1>("Onebody density matrix down", lat,
                                                                       op_vec(1,this->identity_matrix(0)),
                                                                       op_vec(1,this->filling_matrix(0)),
                                                                       ops, false, false) );
        }
        */
        
        return meas;
    }

    tag_type identity_matrix_tag(size_t type) const
    {
        return identity;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill;
    }
    typename SU2U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        typename SU2U1::charge ret(0);
        ret[0] = static_cast<int>(parms["u1_total_charge1"]);
        ret[1] = static_cast<int>(parms["su2_total_charge"]);
        return ret;
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "create")
            return create;
        else if (name == "create_fill")
            return create_fill;
        else if (name == "destroy")
            return destroy;
        else if (name == "destroy_fill")
            return destroy_fill;
        else if (name == "count")
            return count;
        else if (name == "doubly_occ")
            return doubly_occ;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }

    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
private:
    Index<SU2U1> phys;

    Lattice const & lat;
    BaseParameters & parms;

    boost::shared_ptr<TagHandler<Matrix, SU2U1> > tag_handler;
    tag_type create, create_fill, destroy, destroy_fill,
             count, doubly_occ,
             identity, fill;
    

    double get_t (BaseParameters & parms, int i) const
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};

#endif
