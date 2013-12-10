/*****************************************************************************
 *
 * ALPS MPS DMRG Project
 *
 * Copyright (C) 2013 Institute for Theoretical Physics, ETH Zurich
 *               2013-2013 by Alexandr Kosenkov <alex.kosenkov@gmail.com>
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

#ifndef TRACKING_H
#define TRACKING_H

#ifdef AMBIENT_TRACKING
    #define ambient_track(variable) ambient_track_as(variable, #variable)
    #define ambient_track_array(variable, i) ambient_track_as( variable[i], (std::string(#variable) + "[" + std::to_string(i) + "]") )

    template<class Matrix, class SymmGroup> class Boundary;
    template<class Matrix, class SymmGroup> class MPSTensor;
    template<class Matrix, class SymmGroup> class block_matrix;

    template<class Matrix, class SymmGroup> 
    void track(block_matrix<Matrix, SymmGroup>& t, const std::string& label){
        for(int i = 0; i < t.n_blocks(); ++i) 
        ambient::track(t[i], label);
        t.label = label;
    }

    template<class Matrix, class SymmGroup> 
    void track(MPSTensor<Matrix, SymmGroup>& t, const std::string& label){
        track(t.data(), label);
    }

    template<typename T>
    void ambient_track_as(T& object, const std::string& label){
        track(object, label);
    }
#else
    #define ambient_track(variable)
    #define ambient_track_array(variable, i)
    #define ambient_track_as(variable, label)
#endif

#endif
