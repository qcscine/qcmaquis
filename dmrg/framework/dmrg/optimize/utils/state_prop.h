/*****************************************************************************
 *
 * ALPS Project: Algorithms and Libraries for Physics Simulations
 *
 * ALPS Libraries
 *
 * Copyright (C) 2017 by Alberto Baiardi <alberto.baiardi@sns.it>
 *
 * This software is part of the ALPS libraries, published under the ALPS
 * Library License; you can use, redistribute it and/or modify it under
 * the terms of the license, either version 1 or (at your option) any later
 * version.
 *
 * You should have received a copy of the ALPS Library License along with
 * the ALPS Libraries; see the file LICENSE.txt. If not, the license is also
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

#ifndef STATE_PROP_H
#define STATE_PROP_H

#include <dmrg/optimize/ietl_lanczos_solver.h>
#include <ietl/traits.h>

namespace ietl {
    // +--------------------+
    //  STATE_PROP STRUCTURE
    // +--------------------+
    // The state_prop structure contains all the quantites identifying an eigenpair of the local
    // eigenvalue problem
    template<class VS>
    struct state_prop {
        // Types definition
        typedef typename vectorspace_traits<VS>::real_type       real_type;
        typedef typename vectorspace_traits<VS>::scalar_type     scalar_type;
        typedef typename vectorspace_traits<VS>::vector_type     vector_type;
        // Attributes
        vector_type u_, uA_;
        double property_ ;
        real_type theta_ ;
        bool has_property_ ;
        // +------------+
        //  CONSTRUCTORS
        // +------------+
        state_prop() = default ;
        state_prop(vector_type &u, vector_type &uA, real_type &theta) :
                u_(u),
                uA_(uA),
                theta_(theta),
                has_property_(false),
                property_(0.)
        {} ;
        // Copy constructor
        state_prop(state_prop const &rhs) :
            u_(rhs.u_),
            uA_(rhs.uA_),
            theta_(rhs.theta_),
            has_property_(rhs.has_property_),
            property_(rhs.property_)
        {} ;
        // Assigment operator
        state_prop& operator=(state_prop const& rhs)
        {
            if (this != &rhs) {
                u_            = rhs.u_;
                uA_           = rhs.uA_;
                theta_        = rhs.theta_;
                has_property_ = rhs.has_property_ ;
                property_     = rhs.property_ ;
            }
            return *this ;
        }
        // +-------+
        //  METHODS
        // +-------+
        void add_property(double const& prop_toadd)
        {
            if (this->has_property_) {
                throw std::runtime_error("Property already present");
            } else {
                has_property_ = true ;
                property_ = prop_toadd ;
            }
        }
    } ;
}

#endif
