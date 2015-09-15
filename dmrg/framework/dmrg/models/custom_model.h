/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2015 by Michele Dolfi <dolfim@phys.ethz.ch>
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

#ifndef CUSTOM_MODEL_H
#define CUSTOM_MODEL_H

#include "dmrg/models/model.h"


/// Custom Model
/// ------------
/// Model that allows site and bond terms to be added dynamically.
/// DO NOT use this class directly, but use the PIMPL class below, which
/// allows easy conversion to the standard Model class used in other parts
/// of the framework.
template <class Matrix, class SymmGroup>
class custom_model_impl : public model_impl<Matrix, SymmGroup>
{
    typedef model_impl<Matrix, SymmGroup> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::measurements_type measurements_type;
    
public:
    typedef typename base::size_t size_t;
    typedef typename Matrix::value_type value_type;
    typedef typename base::op_t op_t;
    
    custom_model_impl (Index<SymmGroup> const& phys_)
    : tag_handler(new table_type())
    , phys(phys_)
    {
        op_t ident_op = identity_matrix<op_t>(phys);
        ident = tag_handler->register_op(ident_op, tag_detail::bosonic);
    }
    
    void add_siteterm(op_t const& op, size_t i, value_type const& coeff=1.)
    {
        tag_type op_tag = tag_handler->register_op(op, tag_detail::bosonic);
        
        term_descriptor term;
        term.coeff = coeff;
        term.push_back( boost::make_tuple(i, op_tag) );
        this->terms_.push_back(term);
    }
    void add_bondterm(op_t const& op_left, size_t i, op_t const& op_right, size_t j, value_type const& coeff=1.)
    {
        tag_type tag_left  = tag_handler->register_op(op_left, tag_detail::bosonic);
        tag_type tag_right = tag_handler->register_op(op_right, tag_detail::bosonic);
        
        term_descriptor term;
        term.coeff = coeff;
        term.push_back( boost::make_tuple(i, tag_left) );
        term.push_back( boost::make_tuple(j, tag_right) );
        this->terms_.push_back(term);
    }
    
    
    /*** Functions needed my model_impl ***/
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<SymmGroup> const& phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return identity_matrix_tag(type);
    }
    typename SymmGroup::charge total_quantum_numbers(BaseParameters & parms) const
    {
        throw std::runtime_error("Model cannot parse total_quantum_numbers from parameters.");
        return SymmGroup::IdentityCharge;
    }
    measurements_type measurements () const
    {
        return measurements_type();
    }
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        throw std::runtime_error("No operators available in this model.");
        return 0;
    }
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
private:
    boost::shared_ptr<TagHandler<Matrix, SymmGroup> > tag_handler;
    Index<SymmGroup> phys;
    tag_type ident;
};


/// pimpl for CustomModel
template <class Matrix, class SymmGroup>
class CustomModel {
    typedef custom_model_impl<Matrix, SymmGroup> impl_type;
    typedef boost::shared_ptr<impl_type> impl_ptr;
public:
    typedef typename impl_type::size_t size_t;
    typedef typename impl_type::value_type value_type;
    typedef typename impl_type::op_t op_t;
    
    CustomModel(Index<SymmGroup> const& phys)
    : impl_(new impl_type(phys))
    { }
    
    Model<Matrix, SymmGroup> make_model() const
    {
        return Model<Matrix, SymmGroup>(boost::static_pointer_cast<model_impl<Matrix, SymmGroup> >(impl_));
    }
    
    void add_siteterm(op_t const& op, size_t i, value_type const& coeff=1.)
    {
        return impl_->add_siteterm(op, i, coeff);
    }
    void add_bondterm(op_t const& op_left, size_t i, op_t const& op_right, size_t j, value_type const& coeff=1.)
    {
        return impl_->add_bondterm(op_left, i, op_right, j, coeff);
    }
    
private:
    impl_ptr impl_;
};



#endif
