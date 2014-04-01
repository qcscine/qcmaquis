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

#ifndef MODELS_CODED_U1_H
#define MODELS_CODED_U1_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"

/* ****************** HEISENBERG */
template<class Matrix>
class Heisenberg : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

public:   
    Heisenberg (const Lattice& lat, double Jxy, double Jz)
    : tag_handler(new table_type())
    {
        phys.insert(std::make_pair(1, 1));
        phys.insert(std::make_pair(-1, 1));

        op_t ident_op, splus_op, sminus_op, sz_op;
   
        ident_op.insert_block(Matrix(1, 1, 1), -1, -1);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        splus_op.insert_block(Matrix(1, 1, 1), -1, 1);
        
        sminus_op.insert_block(Matrix(1, 1, 1), 1, -1);
        
        sz_op.insert_block(Matrix(1, 1, 0.5), 1, 1);
        sz_op.insert_block(Matrix(1, 1, -0.5), -1, -1);
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(splus,   tag_detail::bosonic)
        REGISTER(sminus,  tag_detail::bosonic)
        REGISTER(sz,      tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    term_descriptor term;
                    term.coeff = Jz;
                    term.push_back( boost::make_tuple(p, sz) );
                    term.push_back( boost::make_tuple(neighs[n], sz) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.coeff = Jxy/2;
                    term.push_back( boost::make_tuple(p, splus) );
                    term.push_back( boost::make_tuple(neighs[n], sminus) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.coeff = Jxy/2;
                    term.push_back( boost::make_tuple(p, sminus) );
                    term.push_back( boost::make_tuple(neighs[n], splus) );
                    this->terms_.push_back(term);
                }
            }
        }
        
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        return;
    }
    
    Index<U1> const& phys_dim(size_t type) const
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
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["u1_total_charge"]);
    }

    measurements_type measurements () const
    {
        return measurements_type();
    }

    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "splus")
            return splus;
        else if (name == "sminus")
            return sminus;
        else if (name == "sz")
            return sz;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }

private:
    Index<U1> phys;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, splus, sminus, sz;
};

/* ****************** HARD CORE BOSONS */
template<class Matrix>
class HCB : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
public:   
    HCB (const Lattice& lat, double t=1)
    : tag_handler(new table_type())
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));

        op_t ident_op;
        op_t create_op, destroy_op, count_op;
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        create_op.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy_op.insert_block(Matrix(1, 1, 1), 1, 0);
        
        count_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(create,  tag_detail::bosonic)
        REGISTER(destroy, tag_detail::bosonic)
        REGISTER(count,   tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, create) );
                    term.push_back( boost::make_tuple(neighs[n], destroy) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, destroy) );
                    term.push_back( boost::make_tuple(neighs[n], create) );
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
    
    Index<U1> const& phys_dim(size_t type) const
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
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["u1_total_charge"]);
    }
    
    measurements_type measurements () const
    {
        return measurements_type();
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n")
            return count;
        else if (name == "bdag")
            return create;
        else if (name == "b")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }

    
private:
    Index<U1> phys;
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, create, destroy, count;
};

/* ****************** BOSE-HUBBARD */
template<class Matrix>
class BoseHubbard : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
    typedef typename Matrix::value_type value_type;
public:
    BoseHubbard (const Lattice& lat_, BaseParameters & model_)
    : lat(lat_)
    , model(model_)
    , tag_handler(new table_type())
    {
        int Nmax = model["Nmax"];
        double t = model["t"];
        double U = model["U"];
        double V = model["V"];
        
        op_t ident_op;
        op_t create_op, destroy_op, count_op, interaction_op;

        phys.insert(std::make_pair(0, 1));
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        
        for (int n=1; n<=Nmax; ++n)
        {
            phys.insert(std::make_pair(n, 1));
            
            ident_op.insert_block(Matrix(1, 1, 1), n, n);
            
            count_op.insert_block(Matrix(1, 1, n), n, n);
            if ((n*n-n) != 0)
                interaction_op.insert_block(Matrix(1, 1, n*n-n), n, n);
            
            
            create_op.insert_block(Matrix(1, 1, std::sqrt(value_type(n))), n-1, n);
            destroy_op.insert_block(Matrix(1, 1, std::sqrt(value_type(n))), n, n-1);
        }
        
        
        /**********************************************************************/
        /*** Create operator tag table ****************************************/
        /**********************************************************************/
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,       tag_detail::bosonic)
        REGISTER(create,      tag_detail::bosonic)
        REGISTER(destroy,     tag_detail::bosonic)
        REGISTER(count,       tag_detail::bosonic)
        REGISTER(interaction, tag_detail::bosonic)
        
#undef REGISTER
        /**********************************************************************/

        
        for (int p=0; p<lat.size(); ++p) {
            /* interaction */
            {
                term_descriptor term;
                term.coeff = U/2.;
                term.push_back( boost::make_tuple(p, interaction) );
                this->terms_.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                /* hopping */
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, create) );
                    term.push_back( boost::make_tuple(neighs[n], destroy) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, destroy) );
                    term.push_back( boost::make_tuple(neighs[n], create) );
                    this->terms_.push_back(term);
                }
                /* nearest-neighborn interaction */
                {
                    term_descriptor term;
                    term.coeff = V;
                    term.push_back( boost::make_tuple(p, count) );
                    term.push_back( boost::make_tuple(neighs[n], count) );
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
    
    Index<U1> const& phys_dim(size_t type) const
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
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["u1_total_charge"]);
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n")
            return count;
        else if (name == "bdag")
            return create;
        else if (name == "b")
            return destroy;
        else if (name == "id")
            return ident;
        else if (name == "fill")
            return ident;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    measurements_type measurements () const
    {
        typedef std::vector<block_matrix<Matrix, U1> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;
        
        if (model["ENABLE_MEASURE[Density]"]) {
            std::string name = "Density";
            meas.push_back( new measurements::average<Matrix, U1>(name, lat,
                                                      op_vec(1,this->identity_matrix(0)),
                                                      op_vec(1,this->filling_matrix(0)),
                                                      op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["ENABLE_MEASURE[Local density]"]) {
            std::string name = "Local density";
            meas.push_back( new measurements::local<Matrix, U1>(name, lat,
                                                    op_vec(1,this->identity_matrix(0)),
                                                    op_vec(1,this->filling_matrix(0)),
                                                    op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["ENABLE_MEASURE[Onebody density matrix]"]) {
            std::string name = "Onebody density matrix";
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create)), false) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy)), false) );
            meas.push_back( new measurements::correlations<Matrix, U1>(name, lat,
                                                           op_vec(1,this->identity_matrix(0)),
                                                           op_vec(1,this->filling_matrix(0)),
                                                           ops, true, false) );
        }
        
        return meas;
    }
    
private:
    const Lattice & lat;
    BaseParameters & model;
    Index<U1> phys;
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident, create, destroy, count, interaction;
};

/* ****************** FREE FERMIONS */
template<class Matrix>
class FreeFermions : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
public:
    FreeFermions (const Lattice& lat, double t=1)
    : lattice(lat)
    , tag_handler(new table_type())
    {
        op_t ident_op;
        op_t create_op, destroy_op, sign_op, dens_op;

        create_op.insert_block(Matrix(1, 1, 1), 0, 1);
        destroy_op.insert_block(Matrix(1, 1, 1), 1, 0);
        
        dens_op.insert_block(Matrix(1, 1, 1), 1, 1);
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
        sign_op.insert_block(Matrix(1, 1, 1), 0, 0);
        sign_op.insert_block(Matrix(1, 1, -1), 1, 1);
        
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        
#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,     tag_detail::bosonic)
        REGISTER(create,    tag_detail::fermionic)
        REGISTER(destroy,   tag_detail::fermionic)
        REGISTER(dens,      tag_detail::bosonic)
        REGISTER(sign,      tag_detail::bosonic)
        
#undef REGISTER

        
        for (int p=0; p<lat.size(); ++p) {
            std::vector<int> neighs = lat.forward(p);
            for (int n=0; n<neighs.size(); ++n) {
                {
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, create) );
                    term.push_back( boost::make_tuple(neighs[n], destroy) );
                    this->terms_.push_back(term);
                }
                {
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.coeff = -t;
                    term.push_back( boost::make_tuple(p, destroy) );
                    term.push_back( boost::make_tuple(neighs[n], create) );
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
    
    Index<U1> const& phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return sign;
    }
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["u1_total_charge"]);
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n")
            return dens;
        else if (name == "cdag")
            return create;
        else if (name == "c")
            return destroy;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    measurements_type measurements () const
    {
        typedef std::vector<block_matrix<Matrix, U1> > op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;
        {
            meas.push_back( new measurements::local<Matrix, U1>("Density", lattice,
                                                    op_vec(1,this->identity_matrix(0)),
                                                    op_vec(1,this->filling_matrix(0)),
                                                    op_vec(1,tag_handler->get_op(dens))) );
        }
        {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(dens)), false) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(dens)), false) );
            meas.push_back( new measurements::correlations<Matrix, U1>("DensityCorrelation", lattice,
                                                           op_vec(1,this->identity_matrix(0)),
                                                           op_vec(1,this->filling_matrix(0)),
                                                           ops, true, false) );
        }
        {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create)), true) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy)), true) );
            meas.push_back( new measurements::correlations<Matrix, U1>("OneBodyDM", lattice,
                                                           op_vec(1,this->identity_matrix(0)),
                                                           op_vec(1,this->filling_matrix(0)),
                                                           ops, true, false) );
        }
        return meas;
    }
    
private:
    Lattice lattice;
    Index<U1> phys;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    tag_type ident;
    tag_type create, destroy, sign, dens;
};

/* ****************** FERMI HUBBARD */
// MD: this model is disabled, because the new Model interface only works for
// Hamiltonians with a single filling_matrix between bond terms.
// For moment the user should jsut use the ALPS Model, which provides the same.
// TODO: Two possibilities
// 1) rewrite the Hamiltonian to have different Jordan-Wiger trasform with only
//    one filling matrix
// 2) unroll the filling matrices explicitly in the bond terms
/*
template<class Matrix>
class FermiHubbardU1 : public model_impl<Matrix, U1> base
{
public:
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
 
    FermiHubbardU1(const Lattice& lat, BaseParameters & parms)
    : tag_handler(new table_type())
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 2));
        phys.insert(std::make_pair(2, 1));
        op_t create_up_op, create_down_op, destroy_up_op, destroy_down_op,
            count_up_op, count_down_op, doubly_occ_op,
            sign_up_op, sign_down_op, fill_op, ident_op;
        
        ident_op.insert_block(Matrix(1, 1, 1), 0, 0);
        ident_op.insert_block(Matrix::identity_matrix(2), 1, 1);
        ident_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        {Matrix tmp(1,2,0); tmp(0,0)=1; create_up_op.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; create_up_op.insert_block(tmp, 1, 2);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; create_down_op.insert_block(tmp, 0, 1);}
        {Matrix tmp(2,1,0); tmp(0,0)=1; create_down_op.insert_block(tmp, 1, 2);}
        
        {Matrix tmp(2,1,0); tmp(0,0)=1; destroy_up_op.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,1)=1; destroy_up_op.insert_block(tmp, 2, 1);}
        {Matrix tmp(2,1,0); tmp(1,0)=1; destroy_down_op.insert_block(tmp, 1, 0);}
        {Matrix tmp(1,2,0); tmp(0,0)=1; destroy_down_op.insert_block(tmp, 2, 1);}
        
        {Matrix tmp(2,2,0); tmp(0,0)=1; count_up_op.insert_block(tmp, 1, 1);}
        count_up_op.insert_block(Matrix(1, 1, 1), 2, 2);
        {Matrix tmp(2,2,0); tmp(1,1)=1; count_down_op.insert_block(tmp, 1, 1);}
        count_down_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        doubly_occ_op.insert_block(Matrix(1, 1, 1), 2, 2);
        
        sign_up_op.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=Matrix::identity_matrix(2); tmp(0,0)=-1; sign_up_op.insert_block(tmp, 1, 1);}
        sign_up_op.insert_block(Matrix(1, 1, -1), 2, 2);
        
        sign_down_op.insert_block(Matrix(1, 1, 1), 0, 0);
        {Matrix tmp=Matrix::identity_matrix(2); tmp(1,1)=-1; sign_down_op.insert_block(tmp, 1, 1);}
        sign_down_op.insert_block(Matrix(1, 1, -1), 2, 2);
        
        gemm(sign_up_op, sign_down_op, fill_op);

#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,         tag_detail::bosonic)
        REGISTER(fill,          tag_detail::bosonic)
        REGISTER(create_up,     tag_detail::fermionic)
        REGISTER(create_down,   tag_detail::fermionic)
        REGISTER(destroy_up,    tag_detail::fermionic)
        REGISTER(destroy_down,  tag_detail::fermionic)
        REGISTER(count_up,      tag_detail::bosonic)
        REGISTER(count_down,    tag_detail::bosonic)
        REGISTER(doubly_occ,    tag_detail::bosonic)
        REGISTER(sign_up,       tag_detail::bosonic)
        REGISTER(sign_down,     tag_detail::bosonic)
        
#undef REGISTER
        
        double U = parms["U"];
        op_t tmp;
        for (int p=0; p<lat.size(); ++p) {
            { // U term
                term_descriptor term;
                term.is_fermionic = false;
                term.coeff = U;
                term.push_back( boost::make_tuple(p, doubly_occ) );
                this->terms_.push_back(term);
            }
            
            std::vector<int> neighs = lat.forward(p);
            for (std::vector<int>::iterator hopto = neighs.begin();
                 hopto != neighs.end(); ++hopto)
            {
                double ti = get_t(parms,
                                  lat.get_prop<int>("type", p, *hopto));
                { // t*cdag_up*c_up
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.fill_operator = sign_up;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(sign_up, create_up);
                    term.coeff = -ti * prod.second;
                    term.push_back( boost::make_tuple(p, prod.first) );
                    term.push_back( boost::make_tuple(*hopto, destroy_up) );
                    this->terms_.push_back(term);
                }
                { // t*c_up*cdag_up
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.fill_operator = sign_up;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(destroy_up, sign_up);
                    term.coeff = -ti * prod.second;
                    term.push_back( boost::make_tuple(p, prod.first) );
                    term.push_back( boost::make_tuple(*hopto, create_up) );
                    this->terms_.push_back(term);
                }
                { // t*cdag_down*c_down
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.fill_operator = sign_down;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(sign_down, create_down);
                    term.coeff = -ti * prod.second;
                    term.push_back( boost::make_tuple(p, prod.first) );
                    term.push_back( boost::make_tuple(*hopto, destroy_down) );
                    this->terms_.push_back(term);
                }
                { // t*c_down*cdag_down
                    term_descriptor term;
                    term.is_fermionic = true;
                    term.fill_operator = sign_down;

                    // Note inverse notation because of notation in operator.
                    std::pair<tag_type, value_type> prod = tag_handler->get_product_tag(destroy_down, sign_down);
                    term.coeff = -ti * prod.second;
                    term.push_back( boost::make_tuple(p, prod.first) );
                    term.push_back( boost::make_tuple(*hopto, create_down) );
                    this->terms_.push_back(term);
                }
            }
        }
    }
    
    Index<U1> const& phys_dim(size_t type) const
    {
        return phys;
    }
    tag_type identity_matrix_tag(size_t type) const
    {
        return ident;
    }
    tag_type filling_matrix_tag(size_t type) const
    {
        return fill;
    }
    typename U1::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return static_cast<int>(parms["u1_total_charge"]);
    }
    
    tag_type get_operator_tag(std::string const & name, size_t type) const
    {
        if (name == "n_up")
            return count_up;
        else if (name == "n_down")
            return count_down;
        else if (name == "cdag_up")
            return create_up;
        else if (name == "cdag_down")
            return create_down;
        else if (name == "c_up")
            return destroy_up;
        else if (name == "c_down")
            return destroy_down;
        else if (name == "id")
            return ident;
        else if (name == "fill")
            return fill;
        else
            throw std::runtime_error("Operator not valid for this model.");
        return 0;
    }
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
    Measurements<Matrix, U1> measurements () const
    {
        return Measurements<Matrix, U1>();
    }
    
private:
    Index<U1> phys;
    tag_type create_up, create_down, destroy_up, destroy_down, count_up, count_down, doubly_occ,
             sign_up, sign_down, fill, ident;

    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    
    double get_t (BaseParameters & parms, int i)
    {
        std::ostringstream key;
        key << "t" << (i+1);
        return (parms.is_set(key.str())) ? parms[key.str()] : parms["t"];
    }
};
*/

#endif
