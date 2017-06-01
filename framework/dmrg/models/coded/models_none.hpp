/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
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

#ifndef MODELS_CODED_NONE_H
#define MODELS_CODED_NONE_H

#include <sstream>

#include "dmrg/models/model.h"
#include "dmrg/models/measurements.h"
#include "dmrg/utils/BaseParameters.h"
#include "dmrg/models/chem/parse_integrals_vib.h"

/* ****************** WATSON HAMILTONIAN */
template<class Matrix>
class WatsonHamiltonianNone : public model_impl<Matrix, TrivialGroup>
{
    typedef model_impl<Matrix, TrivialGroup> base;

    typedef typename base::table_type table_type;
    typedef typename base::table_ptr table_ptr;
    typedef typename base::tag_type tag_type;
    
    typedef typename base::term_descriptor term_descriptor;
    typedef typename base::terms_type terms_type;
    typedef typename base::op_t op_t;
    typedef typename base::measurements_type measurements_type;
    
    typedef typename base::size_t size_t;
    typedef typename Matrix::value_type value_type;
    
public:
    WatsonHamiltonianNone (const Lattice& lat, BaseParameters & model_)
    : model(model_)
    , lattice(lat)
    , tag_handler(new table_type())
    {
        // Model parameters
        typedef Lattice::pos_t pos_t;
        int Nmax = model["Nmax"];
        int L    = model["L"];
        op_t ident_op;
        op_t create_op, destroy_op, count_op;
        TrivialGroup::charge C = TrivialGroup::IdentityCharge;
        size_t N = Nmax+1;
        phys.insert(std::make_pair(C, N));
        ident_op.insert_block(Matrix::identity_matrix(N), C, C);
        Matrix mcount(N,N), mcreate(N,N), mdestroy(N,N);
        for (int n=1; n<=Nmax; ++n) {
            mcount(n,n) = n;
            mcreate(n-1,n)  = std::sqrt(value_type(n));   // input n-1, output n
            mdestroy(n,n-1) = std::sqrt(value_type(n));   // input n,   output n-1
        }
        count_op.insert_block(mcount, C,C);
        create_op.insert_block(mcreate, C,C);
        destroy_op.insert_block(mdestroy, C,C);
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
        // Extraction of the force constants
        alps::numeric::matrix<Lattice::pos_t> idx;
        std::vector<value_type> matrix_elements;
        boost::tie(idx, matrix_elements) = chem_detail::parse_integrals_vib<value_type,TrivialGroup>(model, lat);
        // Definition of the Watson Hamiltonian
        for (std::size_t m=0; m < matrix_elements.size(); ++m) {
            int i = idx(m, 0);
            int j = idx(m, 1);
            int k = idx(m, 2);
            int l = idx(m, 3);
            int h = idx(m, 4);
            int w = idx(m, 5);
            /* Harmonic contribution */
            if ( w==-1 && h==-1 && k==-1 && l==-1 && i!=-1 && j!=-1 ) {
               //if (i==j && i >= 0 && j>= 0) {
               //    // Diagonal harmonic contribution
		       //    {
               //    term_descriptor term;
               //    term.coeff = matrix_elements[m];
               //    term.push_back( boost::make_tuple(i, count) );
               //    this->terms_.push_back(term);
               //    std::cout << term << std::endl;
		       //    }
               //    //
		       //    {
               //    term_descriptor term;
               //    term.coeff = matrix_elements[m]/2.0;
               //    term.push_back( boost::make_tuple(i, ident) );
               //    this->terms_.push_back(term);
               //    std::cout << term << std::endl;
		       //    }
               //}
               //else {
               // Off-diagonal contribution (for calculations not involving NM)
		       value_type scale ;
		       // Scaling factor for G matrix-related terms
 		       if (i < 0 && j < 0) {
		           i = -i-2;
		           j = -j-2;
		           scale = -1.0;
		       }
	           else {
                   scale = 1.0;
               }
               std::vector<pos_t> positions(2);
               std::vector<tag_type> operators(2);
               positions[0] = i;
               positions[1] = j;
               {
               // ++
               operators[0] = create;
               operators[1] = create;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // +-
               operators[0] = create;
               operators[1] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = scale*matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // -+
               operators[0] = destroy;
               operators[1] = create;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = scale*matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // --
               operators[0] = destroy;
               operators[1] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               //}
            }
            /* Cubic contribution */
            else if ( h==-1 && w==-1 && l==-1 && k!=-1 ) {
               std::cout << "Cubic" << std::endl ;
               std::vector<pos_t> positions(3);
               std::vector<tag_type> operators(3);
               positions[0] = i;
               positions[1] = j;
               positions[2] = k;
               {
               // +++
               operators[0] = create;
               operators[1] = create;
               operators[2] = create;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // ++-
               operators[0] = create ;
               operators[1] = create ;
               operators[2] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // +-+
               operators[0] = create ;
               operators[1] = destroy;
               operators[2] = create ;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // -++
               operators[0] = destroy;
               operators[1] = create ;
               operators[2] = create ;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // +--
               operators[0] = create;
               operators[1] = destroy;
               operators[2] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // -+-
               operators[0] = destroy;
               operators[1] = create ;
               operators[2] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // --+
               operators[0] = destroy;
               operators[1] = destroy;
               operators[2] = create;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
               {
               // ---
               operators[0] = destroy;
               operators[1] = destroy;
               operators[2] = destroy;
               term_descriptor term = arrange_operators(positions, operators, tag_handler);
               term.coeff = matrix_elements[m];
               this->terms_.push_back(term);
               std::cout << term << std::endl;
               }
            }
            /* Quartic contribution */
            // For the sake of simplicity, only one term has been added. 
            // All the possible permutations should be included here
            else if ( i!=-1 && j!=-1 && k!=-1 && l!=-1 && h==-1 && w==-1 ) {
                std::vector<pos_t> positions(4);
                std::vector<tag_type> operators(4);
		        value_type scale ;
                std::cout << "Quartic" << std::endl ;
		        // Scaling factor for Coriolis couplings
 		        if (j < 0 && l < 0) {
		            j = -j-2;
		            l = -l-2;
		            scale = -1.0;
		        }
	        	else {
                    scale = 1.0;
                }
                //
                positions[0] = i;
                positions[1] = j;
                positions[2] = k;
                positions[3] = l;
                {
                // ++++
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++-
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-+
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++--
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +---
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+--
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ----
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = scale*matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
            }
            /* Quintic contribution */
            // For the sake of simplicity, only one term has been added. 
            // All the possible permutations should be included here
            else if ( i!=-1 && j!=-1 && k!=-1 && l!=-1 && h!=-1 && w==-1 ) {
                std::vector<pos_t> positions(L+5);
                std::vector<tag_type> operators(L+5);
                std::cout << "Quintic" << std::endl ;
                positions[0] = i;
                positions[1] = j;
                positions[2] = k;
                positions[3] = l;
                positions[4] = h;
                int JCont = 0 ;
                for (int n=0; n<=L-1; ++n){
                    if (n!=i && n!=j && n!=k && n!=l && n!= h){
                        JCont++ ;
                        positions[4+JCont] = n ;
                    }
                }
                positions.resize(JCont+5) ;
                operators.resize(JCont+5) ;
                for  (int n=1; n<=JCont; ++n){
                    operators[n+4] = ident;
                }
                {
                // +++++
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++-+
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-++
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++--+
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+-+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +---+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++-+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+--+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+-+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ----+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++++-
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++--
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-+-
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++---
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-++-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+--
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--+-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +----
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+++-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++--
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-+-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+---
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --++-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+--
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---+-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -----
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
            }
            /* Sextic contribution */
            // For the sake of simplicity, only one term has been added. 
            // All the possible permutations should be included here
            else if ( i!=-1 && j!=-1 && k!=-1 && l!=-1 && h!=-1 && w!=-1 ) {
                std::vector<pos_t> positions(L+6);
                std::vector<tag_type> operators(L+6);
                std::cout << "Sextic" << std::endl ;
                //
                positions[0] = i;
                positions[1] = j;
                positions[2] = k;
                positions[3] = l;
                positions[4] = h;
                positions[5] = w;
                int JCont = 0 ;
                for (int n=0; n<=L-1; ++n){
                    if (n!=i && n!=j && n!=k && n!=l && n!=h && n!=w ){
                        JCont++ ;
                        positions[5+JCont] = n ;
                    }
                }
                positions.resize(JCont+6) ;
                operators.resize(JCont+6) ;
                for  (int n=1; n<=JCont; ++n){
                    operators[n+5] = ident;
                }
                {
                // ++++++
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++-++
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-+++
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++--++
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-++++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+-++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--+++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +---++
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+++++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++-++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-+++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+--++
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --++++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+-++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---+++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ----++
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++++-+
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++--+
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-+-+
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++---+
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-++-+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+--+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--+-+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +----+
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+++-+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++--+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-+-+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+---+
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --++-+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+--+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---+-+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -----+
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = create;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++++-
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++-+-
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-++-
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++--+-
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+++-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+-+-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--++-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +---+-
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++++-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++-+-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-++-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+--+-
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+++-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+-+-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---++-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ----+-
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = create;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++++--
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +++---
                operators[0] = create;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++-+--
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ++----
                operators[0] = create;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-++--
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-+---
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +--+--
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // +-----
                operators[0] = create;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+++--
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -++---
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+-+--
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // -+----
                operators[0] = destroy;
                operators[1] = create;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --++--
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // --+---
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = create;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ---+--
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = create;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
                {
                // ------
                operators[0] = destroy;
                operators[1] = destroy;
                operators[2] = destroy;
                operators[3] = destroy;
                operators[4] = destroy;
                operators[5] = destroy;
                term_descriptor term = arrange_operators(positions, operators, tag_handler);
                term.coeff = matrix_elements[m];
                this->terms_.push_back(term);
                std::cout << term << std::endl;
                }
            }
        }
    }


    term_descriptor
    arrange_operators(std::vector<Lattice::pos_t> const & positions,
                      std::vector<typename OPTable<Matrix, TrivialGroup>::tag_type> const & operators,
                      boost::shared_ptr<TagHandler<Matrix, TrivialGroup> > tag_handler)
    {
        assert(positions.size() == operators.size());
	// Variables definition        
        typedef          Lattice::pos_t                          pos_t;
        typedef typename Matrix::value_type                      value_type;
        typedef typename OPTable<Matrix, TrivialGroup>::tag_type tag_type;
	typedef std::pair<pos_t, tag_type>                       pos_op_t;
        
        std::vector<pos_op_t> pos_ops;
        std::transform(positions.begin(), positions.end(), operators.begin(), std::back_inserter(pos_ops),
                       std::make_pair<pos_t const&, tag_type const&>);
        std::stable_sort(pos_ops.begin(), pos_ops.end(), generate_mpo::compare<pos_op_t>);

        term_descriptor term;
        for (size_t opnr = 0; opnr < pos_ops.size(); )
        {
            tag_type product = pos_ops[opnr].second;
            size_t range_end = opnr + 1;
            while (range_end < pos_ops.size() && pos_ops[range_end].first == pos_ops[opnr].first){
                value_type scale = 1.0;
                boost::tie(product,scale) = tag_handler->get_product_tag(pos_ops[range_end].second,product);
                range_end++;
            }

            term.push_back( boost::make_tuple(pos_ops[opnr].first, product) );

            opnr = range_end;

        } 
        return term;
    }
    
    void update(BaseParameters const& p)
    {
        // TODO: update this->terms_ with the new parameters
        throw std::runtime_error("update() not yet implemented for this model.");
        return;
    }
    
    Index<TrivialGroup> const& phys_dim(size_t type) const
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
    typename TrivialGroup::charge total_quantum_numbers(BaseParameters & parms) const
    {
        return TrivialGroup::IdentityCharge;
    }
    
    
   measurements_type measurements() const
    {
        typedef std::vector<op_t> op_vec;
        typedef std::vector<std::pair<op_vec, bool> > bond_element;
        
        measurements_type meas;
        
        if (model["MEASURE[Density]"]) {
            meas.push_back( new measurements::average<Matrix, TrivialGroup>("Density", lattice,
                                                                  op_vec(1,this->identity_matrix(0)),
                                                                  op_vec(1,this->filling_matrix(0)),
                                                                  op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Local density]"]) {
            meas.push_back( new measurements::local<Matrix, TrivialGroup>("Local density", lattice,
                                                                op_vec(1,this->identity_matrix(0)),
                                                                op_vec(1,this->filling_matrix(0)),
                                                                op_vec(1,tag_handler->get_op(count))) );
        }
        
        if (model["MEASURE[Onebody density matrix]"]) {
            bond_element ops;
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(create)), false) );
            ops.push_back( std::make_pair(op_vec(1,tag_handler->get_op(destroy)), false) );
            meas.push_back( new measurements::correlations<Matrix, TrivialGroup>("Onebody density matrix", lattice,
                                                                                 op_vec(1,this->identity_matrix(0)),
                                                                                 op_vec(1,this->filling_matrix(0)),
                                                                                 ops, true, false) );
        }
        
        return meas;
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
    BaseParameters & model;
    Lattice lattice;
    Index<TrivialGroup> phys;

    table_ptr tag_handler;
    tag_type ident, create, destroy, count;
};

#endif
