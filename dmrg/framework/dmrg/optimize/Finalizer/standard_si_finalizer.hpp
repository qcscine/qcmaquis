/*****************************************************************************
*
* ALPS MPS DMRG Project
*
* Copyright (C) 2014 Institute for Theoretical Physics, ETH Zurich
*               2011-2011 by Bela Bauer <bauerb@phys.ethz.ch>
*               2011-2013    Michele Dolfi <dolfim@phys.ethz.ch>
*               2014-2014    Sebastian Keller <sebkelle@phys.ethz.ch>
*               2019         Leon Freitag <lefreita@ethz.ch>
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
#ifndef STANDARD_SI_FINALIZER_HPP
#define STANDARD_SI_FINALIZER_HPP
template<class Matrix, class VecSpace>
class StandardSIFinalizer : public Finalizer<Matrix, VecSpace> {
public:
    typedef typename ietl::vectorspace_traits<VecSpace>::scalar_type scalar_type;
    typedef typename ietl::number_traits<scalar_type>::magnitude_type magnitude_type;
    typedef typename ietl::vectorspace_traits<VecSpace>::vector_type vector_type;

    magnitude_type compute_energy(size_t i_state, size_t idx)
    {
        magnitude_type res;
        vector_type jnk;
        ietl::mult(*(this->get_Hamiltonian()), (this->get_u(idx)), jnk, i_state);
        res  = std::real(ietl::dot((this->get_u(idx)), jnk)) ;
        res /= std::real(ietl::dot((this->get_u(idx)), (this->get_u(idx))));
        return res;
    }

    magnitude_type theta_converter(magnitude_type theta)
    {
        return std::fabs(1./theta) ;
    }

    vector_type compute_error(size_t i_state, size_t idx)
    {
        vector_type junk ;
        ietl::mult(*(this->get_Hamiltonian()), (this->get_u(idx)), junk, i_state) ;
        junk -= this->compute_energy(i_state, idx) * (this->get_u(idx)) ;
        junk /= ietl::two_norm((this->get_u(idx))) ;
        return junk ;
    }
};
#endif