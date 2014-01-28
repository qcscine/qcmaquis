/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
 *                            Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef MODELS_CODED_2U1_H
#define MODELS_CODED_2U1_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** FERMI HUBBARD */
template<class Matrix>
class FermiHubbardTwoU1 : public model_impl<Matrix, TwoU1>
{
public:
    typedef model_impl<Matrix, TwoU1> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;

    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef typename Matrix::value_type value_type;
    
    FermiHubbardTwoU1(const Lattice& lat_, BaseParameters & parms_)
    : lat(lat_)
    , parms(parms_)
    , tag_handler(new TagHandler<Matrix, TwoU1>())
    {
        TwoU1::charge A(0), B(0), C(0), D(1);
        B[0]=1; C[1]=1;
        phys.insert(std::make_pair(A, 1));
        phys.insert(std::make_pair(B, 1));
        phys.insert(std::make_pair(C, 1));
        phys.insert(std::make_pair(D, 1));
        
        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
             count_up_op, count_down_op, doubly_occ_op,
             fill_op, ident_op;
        
        ident_op.insert_block(Matrix(1, 1, 1), A, A);
        ident_op.insert_block(Matrix(1, 1, 1), B, B);
        ident_op.insert_block(Matrix(1, 1, 1), C, C);
        ident_op.insert_block(Matrix(1, 1, 1), D, D);
        
        create_up_op.insert_block(Matrix(1, 1, 1), A, B);
        create_up_op.insert_block(Matrix(1, 1, 1), C, D);
        create_down_op.insert_block(Matrix(1, 1, 1), A, C);
        create_down_op.insert_block(Matrix(1, 1, 1), B, D);
        
        destroy_up_op.insert_block(Matrix(1, 1, 1), B, A);
        destroy_up_op.insert_block(Matrix(1, 1, 1), D, C);
        destroy_down_op.insert_block(Matrix(1, 1, 1), C, A);
        destroy_down_op.insert_block(Matrix(1, 1, 1), D, B);
        
        count_up_op.insert_block(Matrix(1, 1, 1), B, B);
        count_up_op.insert_block(Matrix(1, 1, 1), D, D);
        count_down_op.insert_block(Matrix(1, 1, 1), C, C);
        count_down_op.insert_block(Matrix(1, 1, 1), D, D);
        
        doubly_occ_op.insert_block(Matrix(1, 1, 1), D, D);
        
        // TODO: use one sign operator only
        fill_op.insert_block(Matrix(1, 1, 1), A, A);
        fill_op.insert_block(Matrix(1, 1, -1), B, B);
        fill_op.insert_block(Matrix(1, 1, -1), C, C);
        fill_op.insert_block(Matrix(1, 1, 1), D, D);

        op_t tmp;

        gemm(fill_op, create_down_op, tmp);
        create_down_op = tmp;
        gemm(destroy_down_op, fill_op, tmp);
        destroy_down_op = tmp; 
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,        tag_detail::bosonic)
        REGISTER(fill,         tag_detail::bosonic)
        REGISTER(create_up,    tag_detail::fermionic)
        REGISTER(create_down,  tag_detail::fermionic)
        REGISTER(destroy_up,   tag_detail::fermionic)
        REGISTER(destroy_down, tag_detail::fermionic)
        REGISTER(count_up,     tag_detail::bosonic)
        REGISTER(count_down,   tag_detail::bosonic)
        REGISTER(doubly_occ,   tag_detail::bosonic)
        
#undef REGISTER
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
                { // t*cdag_up*c_up
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    ptag = tag_handler->get_product_tag(fill, create_up); // Note inverse notation because of notation in operator.
                    term.coeff *= ptag.second;

                    term.push_back( boost::make_tuple(p, ptag.first) );
                    term.push_back( boost::make_tuple(*hopto, destroy_up) );
                    this->terms_.push_back(term);
                }
                { // t*c_up*cdag_up
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    ptag = tag_handler->get_product_tag(fill, destroy_up); // Note inverse notation because of notation in operator.
                    term.coeff *= -ptag.second; // Note minus because of anti-commutation

                    term.push_back( boost::make_tuple(p, ptag.first) );
                    term.push_back( boost::make_tuple(*hopto, create_up) );
                    this->terms_.push_back(term);
                }
                { // t*cdag_down*c_down
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    ptag = tag_handler->get_product_tag(fill, create_down); // Note inverse notation because of notation in operator.
                    term.coeff *= ptag.second;

                    term.push_back( boost::make_tuple(p, ptag.first) );
                    term.push_back( boost::make_tuple(*hopto, destroy_down) );
                    this->terms_.push_back(term);
                }
                { // t*c_down*cdag_down
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    ptag = tag_handler->get_product_tag(fill, destroy_down); // Note inverse notation because of notation in operator.
                    term.coeff *= -ptag.second; // Note minus because of anti-commutation

                    term.push_back( boost::make_tuple(p, ptag.first) );
                    term.push_back( boost::make_tuple(*hopto, create_down) );
                    this->terms_.push_back(term);
                }
            }
        }
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<TwoU1> const & phys_dim(size_t type) const
    {
        return phys;
    }
    
    measurements_type measurements () const
    {
        typedef std::vector<block_matrix<Matrix, TwoU1> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;

        if (parms["ENABLE_MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, TwoU1>("Density up", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_up))) );
        }
        if (parms["ENABLE_MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, TwoU1>("Density down", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_down))) );
        }
        
        if (parms["ENABLE_MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, TwoU1>("Local density up", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_up))) );
        }
        if (parms["ENABLE_MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, TwoU1>("Local density down", lat,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count_down))) );
        }
        
        if (parms["ENABLE_MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create_up)), true) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy_up)), true) );
            meas.push_back( new measurements::correlations<Matrix, TwoU1>("Onebody density matrix up", lat,
                                                                       op_vec(1,this->identity_matrix(0)),
                                                                       op_vec(1,this->filling_matrix(0)),
                                                                       ops, false, false) );
        }
        if (parms["ENABLE_MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create_down)), true) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy_down)), true) );
            meas.push_back( new measurements::correlations<Matrix, TwoU1>("Onebody density matrix down", lat,
                                                                       op_vec(1,this->identity_matrix(0)),
                                                                       op_vec(1,this->filling_matrix(0)),
                                                                       ops, false, false) );
        }
        
        return meas;
    }

    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill;
    }
    typename TwoU1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        typename TwoU1::charge ret(0);
        ret[0] = static_cast<int>(parms["u1_total_charge1"]);
        ret[1] = static_cast<int>(parms["u1_total_charge2"]);
        return ret;
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "create_up")
            return create_up;
        else if (name == "create_down")
            return create_down;
        else if (name == "destroy_up")
            return destroy_up;
        else if (name == "destroy_down")
            return destroy_down;
        else if (name == "count_up")
            return count_up;
        else if (name == "count_down")
            return count_down;
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
    Index<TwoU1> phys;

    Lattice const & lat;
    BaseParameters & parms;

    boost::shared_ptr<TagHandler<Matrix, TwoU1> > tag_handler;
    tag_type create_up, create_down, destroy_up, destroy_down,
             count_up, count_down, doubly_occ,
             ident, fill;
    

    double get_t (BaseParameters & parms, int i) const
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};

#endif
