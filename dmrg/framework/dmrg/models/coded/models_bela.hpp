/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
 *               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
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

#ifndef MODELS_CODED_BELA_H
#define MODELS_CODED_BELA_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/utils/BaseParameters.h"

/* Chiral spins: this one loads adjacency etc from an external file */
template<class Matrix>
class Chiral_ext : public model_impl<Matrix, U1>
{
    typedef model_impl<Matrix, U1> base;
    
    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;

    typedef std::vector<op_t> op_vec;
    
public:
    // This can only work for complex types, that's why there's a callback to do_init with a tag
    Chiral_ext(const Lattice & lat, BaseParameters & parms)
    : tag_handler(new table_type())
    , lat(lat)
    {
        do_init(parms, typename Matrix::value_type());
    }
    
    void do_init(BaseParameters & parms, double)
    {
        phys.insert(std::make_pair(0, 1));
        phys.insert(std::make_pair(1, 1));
        
        maquis::cout << "Wrong guy." << std::endl;
        exit(0);
    }
    
    void ZZ(int p1, int p2)
    {
        std::ostringstream name_oss;
        name_oss << "ZZ " << p1 << " " << p2;
        std::vector< std::pair<int, op_t> > t;
        t.push_back( std::make_pair(p1, sz_op) );
        t.push_back( std::make_pair(p2, sz_op) );
//        measterm.custom_ops.push_back(t);
//        m_terms.push_back(measterm);

        typedef std::vector< std::vector< std::pair<int, op_t> > > craptype;

        meas.push_back( new measurements::custom<Matrix, U1>(name_oss.str(), lat,
                                                             op_vec(1,ident_op),
                                                             op_vec(1,ident_op),
                                                             craptype(1, t) ) );
    }
    
    void do_init(BaseParameters & parms, std::complex<double>)
    {
        if (parms.get<double>("spin") == 1) {
            phys.insert(std::make_pair(-1, 1));
            phys.insert(std::make_pair(0, 1));
            phys.insert(std::make_pair(1, 1));

            ident_op.insert_block(Matrix(1,1,1),-1,-1);
            ident_op.insert_block(Matrix(1,1,1),0,0);
            ident_op.insert_block(Matrix(1,1,1),1,1);

            sz_op.insert_block(Matrix(1, 1, 1), 1, 1);
            sz_op.insert_block(Matrix(1, 1, 0), 0, 0);
            sz_op.insert_block(Matrix(1, 1, -1), -1, -1);

            splus_op.insert_block(Matrix(1, 1, sqrt(2.0)), -1, 0);
            splus_op.insert_block(Matrix(1, 1, sqrt(2.0)), 0, 1);

            sminus_op.insert_block(Matrix(1, 1, sqrt(2.0)), 1, 0);
            sminus_op.insert_block(Matrix(1, 1, sqrt(2.0)), 0, -1);
        } else {
            phys.insert(std::make_pair(-1, 1));
            phys.insert(std::make_pair(1, 1));
            
            ident_op.insert_block(Matrix(1, 1, 1), -1, -1);
            ident_op.insert_block(Matrix(1, 1, 1), 1, 1);
            
            splus_op.insert_block(Matrix(1, 1, 1), -1, 1);
            sminus_op.insert_block(Matrix(1, 1, 1), 1, -1);
            
            sz_op.insert_block(Matrix(1, 1, 0.5), 1, 1);
            sz_op.insert_block(Matrix(1, 1, -0.5), -1, -1);
        }

#define REGISTER(op, kind) op = tag_handler->register_op(op ## _op, kind);
        
        REGISTER(ident,   tag_detail::bosonic)
        REGISTER(splus,   tag_detail::bosonic)
        REGISTER(sminus,  tag_detail::bosonic)
        REGISTER(sz,      tag_detail::bosonic)
        
#undef REGISTER
        
        {
            int p1 = 0, p2;
            for (p2 = 1; p2 < lat.size(); ++p2)
                ZZ(p1,p2);
            
            //for (p1 = 0; p1 < lat.size()-1; ++p1)
            //    ZZ(p1,p1+1);
            
            //p1 = (lat.size()-1)/2;
            //for (p2 = 0; p2 < lat.size(); ++p2)
            //    if (p1 != p2)
            //        ZZ(p1,p2);
           
            int n = (lat.size()-1)/4;
            for (int o = 0; o < n; ++o) {
#define upperleft(i) (2*n-2-2*(i))
#define upperright(i) (2*n+1+2*(i))
#define lowerright(i) (2*n+2+2*(i))
                ZZ(upperleft(o), lowerright(o));
                ZZ(upperleft(o), upperright(o));
                ZZ(upperright(o), lowerright(o));
                if (o < n-1) {
                    ZZ(upperleft(o+1), lowerright(o));
                    ZZ(upperleft(o+1), upperright(o));
                    ZZ(upperright(o+1), lowerright(o));

                    ZZ(upperleft(o), lowerright(o+1));
                    ZZ(upperleft(o), upperright(o+1));
                    ZZ(upperright(o), lowerright(o+1));
                }
#undef upperleft
#undef upperright
#undef lowerright
            }


        }
        
        maquis::cout << "Reading model from " << parms.get<std::string>("mfile") << std::endl;
        std::ifstream modelf( parms.get<std::string>("mfile").c_str() );
        if (!modelf) {
            maquis::cout << "Something wrong with the file. Exiting." << std::endl;
            exit(0);
        }
        
        while (modelf) {
            std::string line;
            getline(modelf, line);
            
            //maquis::cout << line << std::endl;
            
            std::vector<std::string> toks;
            typedef boost::tokenizer<boost::char_separator<char> > tokenizer;
            boost::char_separator<char> sep(" ");
            tokenizer tok(line, sep);
            
            std::copy(tok.begin(), tok.end(), std::back_inserter(toks));
            
            if (toks.size() == 0)
                continue;
            
            if (toks[0] == "t") {
                maquis::cout << "chiral term: ";
                std::copy(toks.begin()+1, toks.end(), maquis::ostream_iterator<std::string>(maquis::cout, " "));
                maquis::cout << std::endl;
                
                int p1 = boost::lexical_cast<int>(toks[1]),
                p2 = boost::lexical_cast<int>(toks[2]),
                p3 = boost::lexical_cast<int>(toks[3]);
                f = std::complex<double>(0, 0.5*boost::lexical_cast<double>(toks[4]));
                
                TERM(p1, p2, p3);
                
            } else if (toks[0] == "h") {
                maquis::cout << "HB term: ";
                std::copy(toks.begin()+1, toks.end(), maquis::ostream_iterator<std::string>(maquis::cout, " "));
                maquis::cout << std::endl;
                
                int p1 = boost::lexical_cast<int>(toks[1]),
                p2 = boost::lexical_cast<int>(toks[2]);
                
                TERM2(p1, p2, boost::lexical_cast<double>(toks[3]));
            } else if (toks[0] == "c") {
                maquis::cout << "CURRENT term: ";
                std::copy(toks.begin()+1, toks.end(), std::ostream_iterator<std::string>(maquis::cout, " "));
                maquis::cout << std::endl;
                
                int p1 = boost::lexical_cast<int>(toks[1]),
                p2 = boost::lexical_cast<int>(toks[2]);
                
                TERM2C(p1, p2, boost::lexical_cast<double>(toks[3]));
            }
        }
        
        int L = lat.size();
        double h = parms.get<double>("h");
        if (h != 0) {
            for (int p = 0; p < L; ++p) {
                maquis::cout << "Field " << h << " on " << p << std::endl;
                term_descriptor term;
                term.coeff = h;
                term.push_back( boost::make_tuple(p, sz) );
                this->terms_.push_back(term);
            }
        }

        meas.push_back( new measurements::local<Matrix, U1>("Sz", lat,
                                                            op_vec(1,ident_op),
                                                            op_vec(1,ident_op),
                                                            op_vec(1,tag_handler->get_op(sz))) );
        
    }
    
    measurements_type measurements () const
    {
        return meas;
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
    
    void update(BaseParameters const& p)
    {
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
    
    table_ptr operators_table() const
    {
        return tag_handler;
    }
    
private:
    void TERM(int p1, int p2, int p3)
    {
#define NNN(p1,op1,p2,op2,p3,op3,val) \
{ term_descriptor term; \
term.coeff = val; \
term.push_back( boost::make_tuple(p1, op1) ); \
term.push_back( boost::make_tuple(p2, op2) ); \
term.push_back( boost::make_tuple(p3, op3) ); \
this->terms_.push_back(term); }
        
//#define NNN_m(p1,op1,p2,op2,p3,op3) \
//{ std::vector< std::pair<int, op_t> > t__; \
//t__.push_back( boost::make_tuple(p1, op1) ); \
//t__.push_back( boost::make_tuple(p2, op2) ); \
//t__.push_back( boost::make_tuple(p3, op3) ); \
//measterm.custom_ops.push_back(t__); }
        
        maquis::cout << "Operator on " << p1 << " " << p2 << " " << p3 << " " << f << std::endl;
        {
            NNN(p1, splus, p2, sminus, p3, sz, f);
            NNN(p1, splus, p2, sz, p3, sminus, -f);
            NNN(p1, sminus, p2, sz, p3, splus, f);
            NNN(p1, sminus, p2, splus, p3, sz, -f);
            NNN(p1, sz, p2, splus, p3, sminus, f);
            NNN(p1, sz, p2, sminus, p3, splus, -f);
        }
//        {
//            mterm_t measterm; \
//            std::ostringstream name_oss; \
//            name_oss << "BE " << p1 << " " << p2 << " " << p3; \
//            measterm.name = name_oss.str(); \
//            measterm.type = mterm_t::Custom; \
//            NNN_m(p1,  f * splus, p2, sminus, p3, sz); \
//            NNN_m(p1, -f * splus, p2, sz, p3, sminus); \
//            NNN_m(p1,  f * sminus, p2, sz, p3, splus); \
//            NNN_m(p1, -f * sminus, p2, splus, p3, sz); \
//            NNN_m(p1,  f * sz, p2, splus, p3, sminus); \
//            NNN_m(p1, -f * sz, p2, sminus, p3, splus); \
//            m_terms.push_back(measterm);
//        }
        
#undef NNN
#undef NNN_m
    }
    
    void TERM2(int p1, int p2, double J)
    {
        
#define NN(p1,op1,p2,op2,val) \
{ term_descriptor term; \
term.coeff = val; \
term.push_back( boost::make_tuple(p1, op1) ); \
term.push_back( boost::make_tuple(p2, op2) ); \
this->terms_.push_back(term); }
        
        maquis::cout << "Operator on " << p1 << " " << p2 << " " << J << std::endl;
        
        NN(p1, splus, p2, sminus, J*0.5);
        NN(p1, sminus, p2, splus, J*0.5);
        NN(p1, sz, p2, sz, J);
        
#undef NN
        
    }
    
    void TERM2C(int p1, int p2, double J)
    {
        
#define NN(p1,op1,p2,op2,val) \
{ term_descriptor term; \
term.coeff = val; \
term.push_back( boost::make_tuple(p1, op1) ); \
term.push_back( boost::make_tuple(p2, op2) ); \
this->terms_.push_back(term); }
        
        maquis::cout << "Current operator on " << p1 << " " << p2 << " " << J << std::endl;
        
        NN(p1, splus, p2, sminus, std::complex<double>(0, J));
        NN(p1, sminus, p2, splus, std::complex<double>(0, -J));
        
#undef NN
        
    }
    
    boost::shared_ptr<TagHandler<Matrix, U1> > tag_handler;
    op_t ident_op, splus_op, sminus_op, sz_op;
    tag_type ident, splus, sminus, sz;
    
    Index<U1> phys;
    Lattice const & lat;
    
    measurements_type meas;

    std::complex<double> f;
};


#endif
