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
class FermiHubbardSU2 : public model_impl<Matrix, TwoU1>
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
    
    FermiHubbardSU2(const Lattice& lat_, BaseParameters & parms_)
    : lat(lat_)
    , parms(parms_)
    , tag_handler(new TagHandler<Matrix, TwoU1>())
    {
        TwoU1::charge A(0), B(0), C(0), D(0);
        A[0] = 2; // 20
        B[0] = 1; B[1] =  1; // 11
        C[0] = 1; C[1] = -1; // 1-1
        // D = 00

        phys.insert(std::make_pair(A,1));
        phys.insert(std::make_pair(B,1));
        phys.insert(std::make_pair(C,1));
        phys.insert(std::make_pair(D,1));

        op_t identity_op;
        identity_op.insert_block(Matrix(1,1,1), A, A);
        identity_op.insert_block(Matrix(1,1,1), B, B);
        identity_op.insert_block(Matrix(1,1,1), C, C);
        identity_op.insert_block(Matrix(1,1,1), D, D);

        op_t fill_ccdag_op;
        fill_ccdag_op.insert_block(Matrix(1,1,1),  A, A);
        fill_ccdag_op.insert_block(Matrix(1,1,1),  D, D);
        fill_ccdag_op.insert_block(Matrix(1,1,-1), B, B); 
        fill_ccdag_op.insert_block(Matrix(1,1,-1), C, C); 
        fill_ccdag_op.insert_block(Matrix(1,1,1),  B, C); 
        fill_ccdag_op.insert_block(Matrix(1,1,1),  C, B); 

        op_t fill_cdagc_op;
        fill_cdagc_op.insert_block(Matrix(1,1,1),  A, A);
        fill_cdagc_op.insert_block(Matrix(1,1,1),  D, D);
        fill_cdagc_op.insert_block(Matrix(1,1,-1), B, B); 
        fill_cdagc_op.insert_block(Matrix(1,1,-1), C, C); 
        fill_cdagc_op.insert_block(Matrix(1,1,-1), B, C); 
        fill_cdagc_op.insert_block(Matrix(1,1,-1), C, B); 

        op_t create_tail_op;
        create_tail_op.insert_block(Matrix(1,1,-2.), B, A);       
        create_tail_op.insert_block(Matrix(1,1,2.), C, A);        
        create_tail_op.insert_block(Matrix(1,1,-sqrt(2.)), D, B); 
        create_tail_op.insert_block(Matrix(1,1,sqrt(2.)), D, C);  

        op_t destroy_tail_op;
        destroy_tail_op.insert_block(Matrix(1,1,1), A, B);        
        destroy_tail_op.insert_block(Matrix(1,1,-1), A, C);       
        destroy_tail_op.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
        destroy_tail_op.insert_block(Matrix(1,1,-sqrt(2.)), C, D);

        op_t create_head_op;
        create_head_op.insert_block(Matrix(1,1,2.), B, A);
        create_head_op.insert_block(Matrix(1,1,2.), C, A);        
        create_head_op.insert_block(Matrix(1,1,sqrt(2.)), D, B);
        create_head_op.insert_block(Matrix(1,1,sqrt(2.)), D, C);  

        op_t destroy_head_op;
        destroy_head_op.insert_block(Matrix(1,1,1), A, B);        
        destroy_head_op.insert_block(Matrix(1,1,1), A, C);        
        destroy_head_op.insert_block(Matrix(1,1,sqrt(2.)), B, D); 
        destroy_head_op.insert_block(Matrix(1,1,sqrt(2.)), C, D); 

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
        
        REGISTER(identity,        tag_detail::bosonic)
        REGISTER(fill_ccdag,   tag_detail::bosonic)
        REGISTER(fill_cdagc,   tag_detail::fermionic)
        REGISTER(create_head,  tag_detail::fermionic)
        REGISTER(create_tail,  tag_detail::fermionic)
        REGISTER(destroy_head, tag_detail::fermionic)
        REGISTER(destroy_tail, tag_detail::bosonic)
        REGISTER(count,        tag_detail::bosonic)
        REGISTER(doubly_occ,   tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

#define PRINT(op) maquis::cout << #op << "\t" << op << std::endl;
        PRINT(identity)
        PRINT(fill_ccdag)  
        PRINT(fill_cdagc)  
        PRINT(create_head) 
        PRINT(create_tail) 
        PRINT(destroy_head)
        PRINT(destroy_tail)
        PRINT(count)       
        PRINT(doubly_occ)  
#undef PRINT
        
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
                { // t*cdag*c
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    for (int fs=0; fs < std::min(*hopto, p); ++fs)
                        term.push_back( boost::make_tuple(fs, identity) );
                    term.push_back( boost::make_tuple(std::min(*hopto,p), create_head) );

                    for (int fs = std::min(*hopto, p)+1; fs < std::max(*hopto, p); ++fs)
                        term.push_back( boost::make_tuple(fs, fill_cdagc) );
                    term.push_back( boost::make_tuple(std::max(*hopto,p), destroy_tail) );

                    for (int fs = std::max(*hopto, p)+1; fs < lat.size(); ++fs)
                        term.push_back( boost::make_tuple(fs, identity) );

                    this->terms_.push_back(term);
                }
                { // t*c*cdag
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -ti;

                    for (int fs=0; fs < std::min(*hopto, p); ++fs)
                        term.push_back( boost::make_tuple(fs, identity) );
                    term.push_back( boost::make_tuple(std::min(*hopto,p), destroy_head) );

                    for (int fs = std::min(*hopto, p)+1; fs < std::max(*hopto, p); ++fs)
                        term.push_back( boost::make_tuple(fs, fill_ccdag) );
                    term.push_back( boost::make_tuple(std::max(*hopto,p), create_tail) );

                    for (int fs = std::max(*hopto, p)+1; fs < lat.size(); ++fs)
                        term.push_back( boost::make_tuple(fs, identity) );

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
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;

        /*
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
        */
        
        return meas;
    }

    tag_type identity_matrix_tag(size_t type) const
    {
        return identity;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill_cdagc;
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
        if (name == "create_head")
            return create_head;
        else if (name == "create_tail")
            return create_tail;
        else if (name == "destroy_head")
            return destroy_head;
        else if (name == "destroy_tail")
            return destroy_tail;
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
    Index<TwoU1> phys;

    Lattice const & lat;
    BaseParameters & parms;

    boost::shared_ptr<TagHandler<Matrix, TwoU1> > tag_handler;
    tag_type create_head, create_tail, destroy_head, destroy_tail,
             count, doubly_occ,
             identity, fill_cdagc, fill_ccdag;
    

    double get_t (BaseParameters & parms, int i) const
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};

#endif
